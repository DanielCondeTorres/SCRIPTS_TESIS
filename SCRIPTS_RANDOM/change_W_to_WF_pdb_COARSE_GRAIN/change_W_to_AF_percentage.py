#!/usr/bin/env python2.7
import subprocess
import optparse
import random
def main() :
    parser = optparse.OptionParser(version='%prog version 1.0')
    parser.add_option('-f', '--infile', help='input pdb file)', action='store')
    options, arguments = parser.parse_args()

#***********************************************************************************************
    #percentage that yo want of AF particles
    percentage=10
    infile =open(str(options.infile),"r")
    ow=[]
    hw1=[]
    hw2=[]
    p=[]; w=[]
    for line in infile:
        if line.strip() :
            fields=line.split()
            if (len(fields)>1) and fields[0]=="ATOM":
                if (fields[2]=="OW"):
                    ow.append(fields[1])
                if (fields[2]=="HW1"):
                    hw1.append(fields[1])
                if (fields[2]=="HW2"):
                    hw2.append(fields[1])
                if (fields[2]=="W"):
                    w.append(fields[1])
    infile.close()
    print('w:',(int(len(w))))
    print ('ow:', (int(len(ow))))
    print ('hw1:', (int(len(hw1))))
    print ('hw2:', (int(len(hw2))))
    print ('lipidos:', (int(len(p))))
   
    number_of_af_particles=int(len(w))*percentage/100
    print('Total to add AF: ',int(number_of_af_particles))
    infile.close()
    #Write new
    count=0;AF=0
    infile =open(str(options.infile),"r")
    infile2=open('new.pdb','w')
    for line in infile:
        line_read=line.split()
        if line_read[0]=='ATOM':
            numero = random.randint(1, 10)

            if line_read[2]=='W' and count<=number_of_af_particles and numero==1:
                line_read[2]='WF';line_read[3]='WF'
                count+=1
                AF+=1
                cfinal="{:6s}{:5s} {:^4s} {:3s}  {:4s}    {:8.6s}{:8.6s}{:8.7s}{:6.4s}{:6.4s} \n".format(line_read[0],line_read[1],line_read[2],line_read[3],line_read[4],line_read[5],line_read[6],line_read[7],line_read[8],line_read[9])
            else:
                cfinal="{:6s}{:5s} {:^4s} {:3s}  {:4s}    {:8.6s}{:8.6s}{:8.7s}{:6.4s}{:6.4s} \n".format(line_read[0],line_read[1],line_read[2],line_read[3],line_read[4],line_read[5],line_read[6],line_read[7],line_read[8],line_read[9])
            infile2.write(cfinal)
        else:
            print(line)
            cfinal=line
            infile2.write(cfinal)

    infile.close()
    infile2.close()
    resultW=subprocess.run(['grep','-wc','W','new.pdb'],capture_output=True)
    resultWF=subprocess.run(['grep','-wc','WF','new.pdb'],capture_output=True)
    print('AF writed: ',int(resultWF.stdout.decode().strip()))
    print('W now: ', int(resultW.stdout.decode().strip()),';; W before: ',int(len(w)))
    print('TOTAL: ',int(resultWF.stdout.decode().strip())+int(resultW.stdout.decode().strip()))
    print('Porcentaje: ', int(resultWF.stdout.decode().strip())/(int(resultWF.stdout.decode().strip())+int(resultW.stdout.decode().strip()))*100)
if __name__=="__main__" :
    main()
