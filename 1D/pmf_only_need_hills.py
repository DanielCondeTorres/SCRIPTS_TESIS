#!/usr/bin/env python3
# python ../../pmf_only_need_hills.py HILLS PMF_FOLDER/
import sys
import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")  # para que funcione en servidores sin pantalla
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def read_fields(hills_file):
    fields = None
    with open(hills_file) as f:
        for line in f:
            if line.startswith("#! FIELDS"):
                fields = line.split()[2:]
                break
    if fields is None:
        raise RuntimeError("No se encontró la línea '#! FIELDS' en el HILLS")
    return fields


def compute_fes_profile(grid, D0, sigma, h, n_hills=None, norm_roi=(0, 7)):
    """
    Devuelve la FES acumulada hasta n_hills.
    Convención: Mínimo en 0 DENTRO de norm_roi.
    """
    if n_hills is None:
        n_hills = len(D0)
    
    F = np.zeros_like(grid)
    # Acumular gaussianas
    for Di, si, hi in zip(D0[:n_hills], sigma[:n_hills], h[:n_hills]):
        F += hi * np.exp(-(grid - Di) ** 2 / (2 * si * si))
    
    Fi = -F # Invertir para FES
    
    # --- Normalización basada en ROI ---
    if norm_roi is not None:
        xmin, xmax = norm_roi
        # Encontrar índices en la grilla que caen dentro del ROI
        roi_indices = np.where((grid >= xmin) & (grid <= xmax))
        
        if roi_indices[0].size > 0:
            # Encontrar el mínimo solo en esa región
            min_in_roi = Fi[roi_indices].min()
            # Restar ese mínimo a *toda* la curva
            Fi -= min_in_roi
        else:
            # Fallback si el ROI está fuera de la grilla (raro)
            Fi -= Fi.min()
    else:
        # Fallback si norm_roi es None (normalización global)
        Fi -= Fi.min()
    # ---
    
    return Fi

# --- FUNCIÓN DE CÁLCULO DE DELTAG (CON INCERTIDUMBRE Y CORRECCIÓN DE ÁREA) ---
def calculate_deltaG_from_areas(
    grid_x,            # D.z coordinates
    fes,               # Normalized FES (min at 0)
    offset_energy,     # Energy of the plateau (plateau_start_y_coord)
    fes_error_estimate, # <<< Incertidumbre de entrada (ej. 2.0 kJ/mol)
    occupied_start,    # Start of the "well" (e.g., 0.0)
    occupied_end,      # End of the "well" (plateau_start_x_coord)
    total_start,       # Start of total integration (e.g., 0.0)
    total_end,         # End of total integration (e.g., 7.0)
    T=298.0,
    R=8.314472/1000.0, # R in kJ/(mol·K)
    out_dir="."
):
    """
    Calcula K, DeltaG, y la incertidumbre de DeltaG (sigma_DeltaG_calc)
    a partir del ratio de áreas de la probabilidad Z(X).
    
    Z_rel(X) = exp( (offset_energy - fes(X)) / RT )
    
    sigma_DeltaG_calc se estima como fes_error_estimate * sqrt(2)
    """
    RT = R * T
    
    # 1. Calcular la probabilidad relativa re-normalizada
    # Usamos np.longdouble para máxima precisión
    fes_ld = fes.astype(np.longdouble)
    offset_energy_ld = np.longdouble(offset_energy)
    RT_ld = np.longdouble(RT)
    Z_rel = np.exp((offset_energy_ld - fes_ld) / RT_ld)
    
    # 2. Definir límites
    a_oc, b_oc = (min(occupied_start, occupied_end), max(occupied_start, occupied_end))
    a_tot, b_tot = (min(total_start, total_end), max(total_start, total_end))

    # 3. Filtrar el dominio al rango total [a_tot, b_tot]
    mask_total = (grid_x >= a_tot) & (grid_x <= b_tot)
    if not np.any(mask_total):
        print("🛑 Dominio total de integración vacío: revisa total_start y total_end.")
        return np.nan, np.nan, np.nan # Devuelve 3 Nans

    X_red = grid_x[mask_total].astype(np.longdouble)
    Z_red = Z_rel[mask_total]
    
    # --- CORRECCIÓN DE ÁREA (calcula áreas por separado) ---

    # 4. Área Ocupada (el "pozo")
    mask_ocupado = (X_red >= a_oc) & (X_red <= b_oc)
    X_ocupado = X_red[mask_ocupado]
    Z_ocupado = Z_red[mask_ocupado]
    area_ocupado = np.trapz(Z_ocupado, X_ocupado) if np.any(mask_ocupado) else 0.0

    # 5. Área Libre (el "plateau")
    mask_libre = (X_red > b_oc) & (X_red <= b_tot)
    X_libre = X_red[mask_libre]
    Z_libre = Z_red[mask_libre]
    area_libre = np.trapz(Z_libre, X_libre) if np.any(mask_libre) else 0.0
    
    # --- FIN DE LA CORRECCIÓN ---

    # 6. Imprimir resultados
    print("-" * 50)
    print("📊 Análisis de DeltaG por Integración de Áreas")
    print(f'Offset Energy (Plateau): {offset_energy:.3f} kJ/mol')
    print(f'Rango Total (Integración): [{a_tot:.3f}, {b_tot:.3f}]')
    print(f'Rango Ocupado (Pozo):     [{a_oc:.3f}, {b_oc:.3f}]')
    # Usar formato 'e' (científico) para manejar números muy grandes o muy pequeños
    print(f'Área Ocupada (Pozo): {area_ocupado:.6e}, Área Libre (Plateau): {area_libre:.6e}')

    # 7. Cálculo de K y DeltaG
    if area_libre <= 1e-12 or area_ocupado <= 1e-12:
        print("🛑 Áreas no positivas o numéricamente nulas → K/ΔG indefinidos.")
        K, DeltaG = np.nan, np.nan
    else:
        K = area_ocupado / area_libre
        DeltaG = -RT_ld * np.log(K)
        print(f'K (Ocupado/Libre): {K:.6e}')
        print(f'DeltaG (calculado): {float(DeltaG):.3f} kJ/mol')
        
    # --- CÁLCULO DE INCERTIDUMBRE ---
    sigma_DeltaG_calc = fes_error_estimate * np.sqrt(2.0)
    print(f'Incertidumbre (σ_ΔG): {sigma_DeltaG_calc:.3f} kJ/mol (basado en σ_FES = {fes_error_estimate:.2f})')
    # ---
    
    print("-" * 50)

    # 8. Visualización (convertir de nuevo a float64 para matplotlib)
    plt.figure(figsize=(8, 5))
    plt.plot(X_red.astype(float), Z_red.astype(float), color='black', lw=2, label=f'Z_rel(X) (Offset={offset_energy:.1f} kJ/mol)')

    # Área ocupada
    plt.fill_between(X_ocupado.astype(float), Z_ocupado.astype(float), alpha=0.5, color='orange', label=f'Área Ocupada (Pozo)\nK = {K:.3e}')

    # Área libre
    plt.fill_between(X_libre.astype(float), Z_libre.astype(float), alpha=0.4, color='skyblue', label='Área Libre (Plateau)')

    # Límites
    plt.axvline(a_oc, color='red', ls='--', label=f'Límite Ocupado Inf. = {a_oc:.3f}')
    plt.axvline(b_oc, color='red', ls='--', label=f'Límite Ocupado Sup. = {b_oc:.3f}')

    plt.xlabel('D.z'); plt.ylabel('Z_rel(X) (Probabilidad relativa re-normalizada)');
    
    # --- Título con error ---
    if np.isnan(DeltaG):
        plt.title(f'Integración de Áreas para DeltaG (Error en cálculo)')
    else:
        plt.title(f'Integración de Áreas para DeltaG\nΔG = {float(DeltaG):.2f} ± {sigma_DeltaG_calc:.2f} kJ/mol')
    # ---
    
    plt.legend(fontsize=8); plt.grid(alpha=0.3); plt.tight_layout()
    
    out_path = os.path.join(out_dir, "Area_Offset_DeltaG.png")
    plt.savefig(out_path, dpi=300)
    plt.close() # Usar close() para script de servidor
    print(f"Guardado: {out_path}")

    return K, DeltaG, sigma_DeltaG_calc # Devuelve 3 valores
# --- FIN DE LA FUNCIÓN DE DELTAG ---

# --- NUEVA FUNCIÓN: BLOCK AVERAGING ---
def estimate_fes_error(grid, all_F_profiles, plot_indices, norm_roi, burn_in_frac=0.5, n_blocks=5):
    """
    Estima la incertidumbre de la FES usando promediado por bloques (Block Averaging).
    
    1. Descarta la fracción 'burn_in_frac' (ej. 0.5 para el 50%) de los perfiles.
    2. Divide los perfiles restantes ("convergidos") en 'n_blocks'.
    3. Calcula la FES promedio para cada bloque.
    4. Normaliza cada FES de bloque (mínimo en ROI = 0).
    5. Calcula la desviación estándar en cada punto de la grilla entre los bloques.
    6. Devuelve el error MÁXIMO (std dev) encontrado dentro del ROI ('plot_indices').
    """
    n_total = len(all_F_profiles)
    start_index = int(n_total * burn_in_frac)
    
    # Perfiles de FES (positivos, crudos) de la parte convergida
    converged_profiles = all_F_profiles[start_index:]
    
    if len(converged_profiles) < n_blocks:
        print(f"Advertencia: No hay suficientes perfiles convergidos ({len(converged_profiles)}) para {n_blocks} bloques.")
        print("Usando error por defecto de 2.0 kJ/mol.")
        return 2.0

    # Dividir los perfiles en bloques
    blocks = np.array_split(converged_profiles, n_blocks)
    
    normalized_block_fes = []
    
    for block in blocks:
        if len(block) == 0:
            continue
            
        # 1. Calcular el promedio del bloque (sigue siendo F positivo)
        avg_F_block = np.mean(block, axis=0)
        
        # 2. Normalizar este promedio (invertir, min en ROI = 0)
        avg_fes_block = -avg_F_block
        
        # Encontrar mínimo en el ROI
        min_in_roi = avg_fes_block[plot_indices].min()
        avg_fes_block_norm = avg_fes_block - min_in_roi
        
        normalized_block_fes.append(avg_fes_block_norm)

    if len(normalized_block_fes) < 2:
        print("Advertencia: No se pudieron crear al menos 2 bloques para calcular std dev.")
        print("Usando error por defecto de 2.0 kJ/mol.")
        return 2.0

    # 3. Apilar los perfiles de bloque normalizados
    stacked_blocks = np.array(normalized_block_fes)
    
    # 4. Calcular la desviación estándar entre los bloques (en cada punto de la grilla)
    std_dev_profile = np.std(stacked_blocks, axis=0)
    
    # 5. Encontrar el error máximo dentro del ROI
    max_error_in_roi = np.max(std_dev_profile[plot_indices])
    
    return max_error_in_roi
# --- FIN DE LA FUNCIÓN BLOCK AVERAGING ---


def make_movie_og(hills_file, out_dir, fields, data, D0, sigma, h):
    """
    Genera un vídeo MP4 de la FES acumulada (solo la curva más reciente).
    (Esta función usa tiempo absoluto, 0 -> final)
    """
    try:
        idx_time = fields.index("time")
    except ValueError:
        print("No hay columna 'time' en el HILLS; no se puede hacer el vídeo.")
        return

    time_ps = data[:, idx_time]
    time_ns = time_ps / 1000.0
    
    total = len(D0)
    ROI_X = (0, 7) 
    Dmin, Dmax = D0.min() - 2 * sigma.max(), D0.max() + 2 * sigma.max()
    grid = np.linspace(Dmin, Dmax, 400)
    
    roi_indices = np.where((grid >= ROI_X[0]) & (grid <= ROI_X[1]))
    if roi_indices[0].size == 0:
        roi_indices = np.where((grid >= grid.min()) & (grid <= grid.max()))

    F_final = compute_fes_profile(grid, D0, sigma, h, n_hills=total, norm_roi=ROI_X)
    global_max = F_final[roi_indices].max() 

    N_FRAMES = min(200, total)
    frame_indices = np.linspace(0, total - 1, N_FRAMES, dtype=int)

    frames_dir = os.path.join(out_dir, "_frames_fes_movie")
    os.makedirs(frames_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(5, 5))
    line, = ax.plot([], [], lw=1)
    
    ax.set_xlim(ROI_X[0], ROI_X[1])
    ax.set_ylim(bottom=0, top=global_max * 1.05)
    
    ax.set_xlabel("D.z")
    ax.set_ylabel("Energy (kJ/mol)")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    title = ax.set_title("")

    print(f"Generando {N_FRAMES} frames en {frames_dir} ...")

    for k, i in enumerate(frame_indices):
        Fi = compute_fes_profile(grid, D0, sigma, h, n_hills=i + 1, norm_roi=ROI_X)
        line.set_data(grid, Fi)
        title.set_text(f"HILLS – hill {i+1}/{total} (time = {time_ns[i]:.1f} ns)")
        frame_path = os.path.join(frames_dir, f"frame_{k:05d}.png")
        fig.savefig(frame_path, dpi=250)

    plt.close(fig)
    
    out_movie = os.path.join(out_dir, "fes_movie.mp4")
    cmd = [
        "ffmpeg", "-y", "-framerate", "24",
        "-i", os.path.join(frames_dir, "frame_%05d.png"),
        "-qscale:v", "1", out_movie,
    ]
    print("Llamando a ffmpeg para crear el vídeo...")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Vídeo creado: {out_movie}")
    except subprocess.CalledProcessError as e:
        print("⚠ ffmpeg ha fallado; se conservan los PNG en", frames_dir)
        print(e.stderr)


def make_movie(hills_file, out_dir, fields, data, D0, sigma, h):
    """
    Genera un vídeo de la FES acumulada (estilo gradiente)
    CON BARRA DE COLOR en ns (tiempo absoluto, 0 -> final)
    """
    try:
        idx_time = fields.index("time")
    except ValueError:
        print("No hay columna 'time' en el HILLS; no se puede hacer el vídeo de gradiente.")
        return

    time_ps = data[:, idx_time]
    time_ns = time_ps / 1000.0
    total = len(D0)
    
    t_min_val = time_ns[0]
    t_max_val = time_ns[-1]
    cbar_label = "Time (ns)"
    
    ROI_X = (0, 7) 
    Dmin, Dmax = D0.min() - 2 * sigma.max(), D0.max() + 2 * sigma.max()
    grid = np.linspace(Dmin, Dmax, 400)
    
    roi_indices = np.where((grid >= ROI_X[0]) & (grid <= ROI_X[1]))
    if roi_indices[0].size == 0:
        roi_indices = np.where((grid >= grid.min()) & (grid <= grid.max()))

    N_FRAMES = min(200, total)
    frame_indices = np.linspace(0, total - 1, N_FRAMES, dtype=int)

    profiles = []
    for i in frame_indices:
        Fi = compute_fes_profile(grid, D0, sigma, h, n_hills=i + 1, norm_roi=ROI_X)
        profiles.append(Fi)
    profiles = np.array(profiles)

    global_max = profiles[:, roi_indices[0]].max()

    frames_dir = os.path.join(out_dir, "_frames_fes_movie_grad")
    os.makedirs(frames_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(5.5, 5), dpi=200)
    cmap = plt.get_cmap("Blues")
    
    norm = plt.Normalize(vmin=t_min_val, vmax=t_max_val)
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    print(f"Generando {N_FRAMES} frames en {frames_dir} ...")

    for k in range(N_FRAMES):
        ax.clear()
        ax.set_xlim(ROI_X[0], ROI_X[1])
        ax.set_ylim(bottom=0, top=global_max * 1.05)
        ax.set_xlabel("D.z")
        ax.set_ylabel("Energy (kJ/mol)")
        i_hill = frame_indices[k]
        # --- CORREGIDO (i_hill en lugar de i+hill) ---
        ax.set_title(f"HILLS – hill {i_hill+1}/{total} (time = {time_ns[i_hill]:.1f} ns)")

        for j in range(k + 1):
            Fi = profiles[j]
            current_hill_index = frame_indices[j]
            current_time = time_ns[current_hill_index]
            color = cmap(norm(current_time))
            alpha_frac = j / (N_FRAMES - 1) 
            alpha = 0.2 + 0.8 * alpha_frac
            rgba = list(color)
            rgba[-1] = alpha
            ax.plot(grid, Fi, color=rgba, linewidth=1)

        cbar = fig.colorbar(sm, ax=ax, label=cbar_label)
        frame_path = os.path.join(frames_dir, f"frame_{k:05d}.png")
        fig.savefig(frame_path, dpi=fig.dpi)
        cbar.remove()

    plt.close(fig)

    out_movie_mp4 = os.path.join(out_dir, "fes_movie_gradient.mp4")
    cmd_mp4 = [
        "ffmpeg", "-y", "-framerate", "15",
        "-i", os.path.join(frames_dir, "frame_%05d.png"),
        "-qscale:v", "1", out_movie_mp4,
    ]
    print("Llamando a ffmpeg para crear el vídeo MP4...")
    try:
        subprocess.run(cmd_mp4, check=True, capture_output=True, text=True)
        print(f"Vídeo MP4 creado: {out_movie_mp4}")
    except subprocess.CalledProcessError as e:
        print("⚠ ffmpeg ha fallado al crear MP4.")
        print(e.stderr)


def main():
    if not (3 <= len(sys.argv) <= 4):
        print("Uso: python hills_fes.py /ruta/al/HILLS /ruta/salida [movie]")
        sys.exit(1)

    hills_file = sys.argv[1]
    out_dir = sys.argv[2]
    make_movie_flag = len(sys.argv) == 4 and sys.argv[3].lower() in ("movie", "--movie")

    if not os.path.isfile(hills_file):
        print(f"ERROR: no existe el archivo HILLS: {hills_file}")
        sys.exit(1)

    os.makedirs(out_dir, exist_ok=True)
    
    ROI_X = (0, 7) # Region Of Interest
    
    # --- Límite inferior para la integración del "pozo" ---
    OCCUPIED_START_X = 0.0 # <-- ¡Elige el inicio del "pozo" (ocupado) aquí!
    
    # --- (Estimación de error AHORA SE CALCULA AUTOMÁTICAMENTE) ---
    # FES_ERROR_ESTIMATE = 2.0 # <-- Eliminado

    # ---- leer HILLS ----
    fields = read_fields(hills_file)
    data = np.loadtxt(hills_file)

    try:
        idx_Dz = fields.index("D.z")
        idx_sigma = fields.index("sigma_D.z")
        idx_h = fields.index("height")
    except ValueError as e:
        raise RuntimeError(f"Faltan columnas necesarias en HILLS: {e}")

    D0 = data[:, idx_Dz]
    sigma = data[:, idx_sigma]
    h = data[:, idx_h]

    # ---- grilla en D.z ----
    Dmin, Dmax = D0.min() - 2 * sigma.max(), D0.max() + 2 * sigma.max()
    grid = np.linspace(Dmin, Dmax, 400)
    
    roi_indices = np.where((grid >= ROI_X[0]) & (grid <= ROI_X[1]))
    if roi_indices[0].size == 0:
        plot_indices = np.arange(len(grid))
    else:
        # Asegurarnos que plot_indices es un array simple
        plot_indices = roi_indices[0] 

    # ---- FES acumulada y snapshots ----
    
    # --- Lógica para la "cuenta atrás" de 100 ns ---
    total = len(D0)
    time_ps = None # Inicializar
    try:
        idx_time = fields.index("time")
        time_ps = data[:, idx_time]
        time_ns = time_ps / 1000.0 
        
        time_total_ns = time_ns[-1]
        time_start_window_ns = time_total_ns - 100.0 # 100 ns
        
        start_index_100ns = np.searchsorted(time_ns, time_start_window_ns)
        start_index_100ns = max(0, start_index_100ns)
        
        t_min_cbar = -100.0  # El más antiguo (claro)
        t_max_cbar = 0.0     # El más reciente (oscuro)
        cbar_label = "Time before end (ns)"
        
        N_HILLS_100NS = total - start_index_100ns
        print(f"Detectados {N_HILLS_100NS} hills en los últimos 100 ns.")
        START_IDX_SNAPSHOTS = start_index_100ns
        TITLE_FIG_2 = "HILLS - Convergence (last 100 ns)"
        OUT_FIG_2 = "fes_convergence_100ns.png"

    except (ValueError, IndexError):
        # Fallback (cuenta atrás de Hill Index)
        print("No se encontró 'time', usando N_LAST=20000 hills por defecto.")
        N_LAST_FALLBACK = 20000
        if N_LAST_FALLBACK > total:
            N_LAST_FALLBACK = total
            
        START_IDX_SNAPSHOTS = total - N_LAST_FALLBACK
        
        t_min_cbar = -N_LAST_FALLBACK # El más antiguo (claro)
        t_max_cbar = 0                 # El más reciente (oscuro)
        cbar_label = "Hills before end"
        
        TITLE_FIG_2 = f"HILLS - Convergence (last {N_LAST_FALLBACK} hills)"
        OUT_FIG_2 = f"fes_convergence_last_{N_LAST_FALLBACK}.png"
    # ---

    F = np.zeros_like(grid) 
    snapshots = [] # Para Fig 2 (últimos 100 ns)
    
    # --- MODIFICADO: Guardar TODOS los perfiles F para el Block Averaging ---
    all_F_profiles = [] # Lista para guardar todos los perfiles F (positivos)
    # ---

    for i, (Di, si, hi) in enumerate(zip(D0, sigma, h)):
        F += hi * np.exp(-(grid - Di) ** 2 / (2 * si * si))
        
        # Guardar todos los perfiles F (positivos)
        all_F_profiles.append(F.copy()) 
        
        # Guardar solo los snapshots de los últimos 100 ns para Fig 2
        if i >= START_IDX_SNAPSHOTS:
            snapshots.append(F.copy())

    # --- MODIFICADO: Estimar el error ANTES de hacer las figuras ---
    print("\n--- Estimando error de FES por Block Averaging ---")
    FES_ERROR_ESTIMATE = estimate_fes_error(
        grid=grid,
        all_F_profiles=all_F_profiles,
        plot_indices=plot_indices, # Pasar los índices del ROI
        norm_roi=ROI_X
    )
    print(f"Error de FES estimado (max std dev in ROI): {FES_ERROR_ESTIMATE:.2f} kJ/mol")
    # ---
            
    # ---- Análisis de la FES Final (para líneas) ----
    # Usar el último perfil de la lista
    F_final_raw = -all_F_profiles[-1] 
    min_val_in_roi = F_final_raw[plot_indices].min()
    F_final_norm = F_final_raw - min_val_in_roi

    # --- Lógica de derivadas para encontrar líneas ---
    
    # 1. Encontrar el Mínimo Global (dentro del ROI)
    min_idx_local = F_final_norm[plot_indices].argmin()
    min_idx_global = plot_indices[min_idx_local]
    min_x_coord = grid[min_idx_global]
    min_y_coord = F_final_norm[min_idx_global]

    # 2. Calcular derivada
    grid_roi = grid[plot_indices]
    fes_roi = F_final_norm[plot_indices]
    dF_roi = np.gradient(fes_roi, grid_roi[1] - grid_roi[0])
    
    plateau_start_x_coord = None
    plateau_start_y_coord = None

    # 3. Buscar el punto de inflexión (pendiente máxima) DESPUÉS del mínimo
    search_range_after_min = np.arange(min_idx_local, len(dF_roi))
    
    if search_range_after_min.size > 0:
        dF_after_min = dF_roi[search_range_after_min]
        positive_slopes_after_min_idx = np.where(dF_after_min > 0)[0]

        if positive_slopes_after_min_idx.size > 0:
            inflection_idx_local_in_range = positive_slopes_after_min_idx[dF_after_min[positive_slopes_after_min_idx].argmax()]
            inflection_idx_local_roi = search_range_after_min[inflection_idx_local_in_range]
            max_slope = dF_roi[inflection_idx_local_roi]
            
            # 4. Buscar el inicio del plateau DESPUÉS del p. de inflexión
            threshold_slope = max_slope * 0.05 
            search_range_for_plateau = np.arange(inflection_idx_local_roi, len(dF_roi))
            
            if search_range_for_plateau.size > 0:
                dF_in_plateau_range = dF_roi[search_range_for_plateau]
                plateau_start_local_indices = np.where(dF_in_plateau_range < threshold_slope)[0]
                
                if plateau_start_local_indices.size > 0:
                    plateau_start_idx_local_in_range = plateau_start_local_indices[0]
                    plateau_start_idx_local_roi = search_range_for_plateau[plateau_start_idx_local_in_range]
                    plateau_start_global_idx = plot_indices[plateau_start_idx_local_roi]
                    
                    plateau_start_x_coord = grid[plateau_start_global_idx]
                    plateau_start_y_coord = F_final_norm[plateau_start_global_idx]

    # Imprimir resultados
    print(f"Mínimo encontrado en D.z = {min_x_coord:.2f} (E = {min_y_coord:.1f})")
    if plateau_start_x_coord is not None:
        print(f"Inicio del plateau (pendiente < 5% max) en D.z = {plateau_start_x_coord:.2f} (E = {plateau_start_y_coord:.1f})")
    else:
        print("No se pudo detectar el inicio del plateau con la lógica de derivadas.")
    # ---

    # ---- figura 1: FES final ----
    
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(grid, F_final_norm)
    ax.set_xlabel("D.z")
    ax.set_ylabel("Energy (kJ/mol)")
    ax.set_title("HILLS - Final FES")
    
    ax.set_xlim(ROI_X[0], ROI_X[1])
    plot_max = fes_roi.max() 
    ax.set_ylim(bottom=0, top=plot_max * 1.05) 
    
    # --- Añadir líneas y texto ---
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()
    
    # Texto de D.z (línea roja) arriba
    ax.axvline(x=min_x_coord, color='r', linestyle='--', linewidth=1)
    ax.text(min_x_coord, ymax * 0.95, f'D.z = {min_x_coord:.2f}', color='r', ha='center', va='top', fontsize=9)

    if plateau_start_y_coord is not None:
        # Texto de Energía (línea negra) a la izquierda
        ax.axhline(y=plateau_start_y_coord, color='k', linestyle='--', linewidth=1)
        ax.text(xmin + (xmax - xmin) * 0.02, plateau_start_y_coord + (ymax - ymin) * 0.01, # Ligeramente encima
                f'Plateau E = {plateau_start_y_coord:.1f} kJ/mol', color='k', ha='left', va='bottom', fontsize=9)

        # Texto de D.z (línea naranja) arriba
        ax.axvline(x=plateau_start_x_coord, color='orange', linestyle='--', linewidth=1)
        ax.text(plateau_start_x_coord, ymax * 0.90, f'Plateau Start D.z = {plateau_start_x_coord:.2f}', color='orange', ha='center', va='top', fontsize=9)
    # ---
    
    fig.tight_layout()
    out_fes = os.path.join(out_dir, "fes.png")
    fig.savefig(out_fes, dpi=600)
    plt.close(fig)
    print(f"Guardado: {out_fes}")

    # ---- figura 2: "Stride" de los últimos N (100 ns) ----
    
    N_IN_SNAPSHOT = len(snapshots)
    
    if N_IN_SNAPSHOT == 0:
        print("No hay snapshots para la figura 2, saltando.")
    else:
        N_STRIDES = 100 
        if N_IN_SNAPSHOT < N_STRIDES:
            N_STRIDES = N_IN_SNAPSHOT
            
        indices_to_plot = np.linspace(0, N_IN_SNAPSHOT - 1, N_STRIDES, dtype=int)
        
        cmap = cm.get_cmap("Blues") 
        fig, ax = plt.subplots(figsize=(5.5, 5))

        all_snaps_norm = []
        for snap_idx in indices_to_plot:
            Fi_raw = -snapshots[snap_idx]
            min_in_roi = Fi_raw[plot_indices].min()
            Fi_norm = Fi_raw - min_in_roi
            all_snaps_norm.append(Fi_norm)
        
        if not all_snaps_norm:
             global_max_fig2 = 1.0 
        else:
            global_max_fig2 = max(f[plot_indices].max() for f in all_snaps_norm)

        for i, Fi_norm in enumerate(all_snaps_norm):
            frac = i / (N_STRIDES - 1) 
            color = cmap(frac)
            alpha = 0.2 + 0.8 * frac
            rgba = list(color)
            rgba[-1] = alpha
            ax.plot(grid, Fi_norm, color=rgba, linewidth=1)

        ax.set_xlabel("D.z")
        ax.set_ylabel("Energy (kJ/mol)")
        ax.set_title(TITLE_FIG_2)
        ax.set_xlim(ROI_X[0], ROI_X[1])
        ax.set_ylim(bottom=0, top=global_max_fig2 * 1.05)
        
        # --- Añadir líneas y texto ---
        ymin, ymax = ax.get_ylim()
        xmin, xmax = ax.get_xlim()

        # Texto de D.z (línea roja) arriba
        ax.axvline(x=min_x_coord, color='r', linestyle='--', linewidth=1.5)
        ax.text(min_x_coord, ymax * 0.95, f'D.z = {min_x_coord:.2f}', color='r', ha='center', va='top', fontsize=9)

        if plateau_start_y_coord is not None:
            # Texto de Energía (línea negra) a la izquierda
            ax.axhline(y=plateau_start_y_coord, color='k', linestyle='--', linewidth=1.5)
            ax.text(xmin + (xmax - xmin) * 0.02, plateau_start_y_coord + (ymax - ymin) * 0.01,
                    f'Plateau E = {plateau_start_y_coord:.1f} kJ/mol', color='k', ha='left', va='bottom', fontsize=9)

            # Texto de D.z (línea naranja) arriba
            ax.axvline(x=plateau_start_x_coord, color='orange', linestyle='--', linewidth=1.5)
            ax.text(plateau_start_x_coord, ymax * 0.90, f'Plateau Start D.z = {plateau_start_x_coord:.2f}', color='orange', ha='center', va='top', fontsize=9)
        # ---
        
        # --- Añadir el colorbar (cuenta atrás) ---
        norm = plt.Normalize(vmin=t_min_cbar, vmax=t_max_cbar)
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label=cbar_label)
        # ---
        
        fig.tight_layout() 
        out_last = os.path.join(out_dir, OUT_FIG_2)
        fig.savefig(out_last, dpi=600)
        plt.close(fig) 
        print(f"Guardado: {out_last}")


    # ---- vídeo opcional ----
    if make_movie_flag and time_ps is not None: # Solo hacer vídeo si tenemos 'time'
        print("\n--- Creando vídeo (estilo gradiente con colorbar) ---")
        make_movie(hills_file, out_dir, fields, data, D0, sigma, h)
        print("\n--- Creando vídeo (estilo simple) ---")
        make_movie_og(hills_file, out_dir, fields, data, D0, sigma, h)
    elif make_movie_flag:
        print("\nNo se puede crear el vídeo, falta la columna 'time' en el fichero HILLS.")


    # --- Llamada a la función de DeltaG (con error) ---
    if plateau_start_x_coord is not None and plateau_start_y_coord is not None:
        print("\n--- Calculando DeltaG por integración de áreas ---")
        K_calc, DeltaG_calc, Sigma_calc = calculate_deltaG_from_areas(
            grid_x=grid,
            fes=F_final_norm,
            offset_energy=plateau_start_y_coord,
            fes_error_estimate=FES_ERROR_ESTIMATE, # <-- Pasa el error calculado
            occupied_start=OCCUPIED_START_X,     # Límite inferior del pozo (0.0 por defecto)
            occupied_end=plateau_start_x_coord,  # Límite superior del pozo (inicio del plateau)
            total_start=OCCUPIED_START_X,        # Límite inferior total (como pediste)
            total_end=ROI_X[1],                  # Límite superior total (7.0)
            out_dir=out_dir
        )
        
        if not np.isnan(DeltaG_calc):
            print(f"\nResultado Final: ΔG = {float(DeltaG_calc):.2f} ± {Sigma_calc:.2f} kJ/mol")
        
    else:
        print("\nNo se pudo calcular DeltaG por integración, faltan datos del plateau.")
    # ---


if __name__ == "__main__":
    main()
