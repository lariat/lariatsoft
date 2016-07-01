import os

#txtfile = open("Short.txt")
txtfile = open("LArIAT_Lifetime_Run1_withErrors.txt")


eraNumb = -1
lowerBoundRun = []
upperBoundRun = []

run_num  = []
lifetime = []
for line in txtfile:
    column = line.split()
    run_num.append(int(column[0]))
    lifetime.append(float(column[1]))

if not len(run_num) == len(lifetime):
    print "ufffffffffff"

for i in xrange(len(run_num)):
    if i == 0: 
        lowerBoundRun.append(0) 
        upperBoundRun.append(run_num[i+1] -1)
        
    elif i == len(run_num) - 1:
        lowerBoundRun.append(run_num[i]) 
        upperBoundRun.append(99999)
    else:
        lowerBoundRun.append(run_num[i]) 
        upperBoundRun.append(run_num[i+1] -1)
             

for i in xrange(len(run_num)):
    EraNum = str(i)
    ############################ Creating New Directory ##############################
    dir_name = "Era"+ EraNum + "_runRange_"+ str(lowerBoundRun[i]) + "_" + str(upperBoundRun[i]) + "_lifetime_" + str(lifetime[i])
    make_dir = "mkdir " + dir_name
    os.system(make_dir)
    
    #################### Creating New fcl and Lifetime Change #######################
    fclName      = dir_name + "/Reco_"+ dir_name  +".fcl"
    cp_fcl       = "cp Reco.fcl " + fclName
    sub_lifetime = "sed -i \'s/LIFETIMETOBECHANGED/"+ str(lifetime[i]) +"/g\' " + fclName
    os.system(cp_fcl)
    os.system(sub_lifetime)

    ############################## Creating New xml #################################

    xmlName      = dir_name + "/ProdReco_"+ dir_name  +".xml"
    cp_xml       = "cp ProdReco.xml " + xmlName
    sub_name     = "sed -i \'s/XXX/Era_"+ EraNum +"/g\' " + xmlName
    os.system(cp_xml)
    os.system(sub_name)

    ############################## Creating A little log file  #################################
    filename = dir_name + "/LittleLog_Era_"+ EraNum +".log"
    target = open(filename, 'w')
    line1 = "file_range: [ " + str(lowerBoundRun[i])  + " , "+ str(upperBoundRun[i])  + " ]"
    line2 = "lifetime: "+ str(lifetime[i])
    target.write(line1)
    target.write("\n")
    target.write(line2)
    target.close()

    ################################ Print out SAM defintions  ##################################
    samDef = open("samDefitionsFile.txt", 'a')
    line3 = "samweb create-definition RecoEra"+ EraNum + "_lovely \"defname: Lovely1 and data_tier digits and run_number > " + str(lowerBoundRun[i]) +" and run_number < "+ str(upperBoundRun[i]) +" \""
    samDef.write(line3)
    samDef.write("\n")
    





