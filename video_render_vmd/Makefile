crear_video:
	ffmpeg -framerate 10 -pattern_type glob -i '*.ppm' -c:v libx264 -crf 1 -movflags +faststart output.mkv
clean:
	rm *ppm #output*

#FORMAT: MPEG-1
#RENDER: INTERNAL TACHYON
#MOVIE SETTINGS: TRAJECTORY
