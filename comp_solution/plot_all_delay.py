import os
import numpy as np
import matplotlib.pyplot as plt
from losoto.h5parm import h5parm

### Loading data

#casa
w, x, y, z = np.loadtxt('delay.txt', delimiter=' ', usecols = (0, 1, 2, 3), unpack = True)

#aips
f, t, a, p1, p2 = np.genfromtxt('sn.txt', usecols = (0, 1, 3, 4, 5), dtype = str, unpack = True)

frequency = []
time = []
antenna = []
pol1 = []
pol2 = []

for i in range(len(t)):
    frequency.append(float(f[i]))
    time.append(float(t[i]))
    antenna.append(a[i])
    pol1.append(float(p1[i]))
    pol2.append(float(p2[i]))


#h5
name_h5="/data020/scratch/lb_bw/test_losoto_delay_calibration/16sec_1SB_2comp_uvcutCS/ILTJ1327_16sec_1SB_2comp_uvcutCS_globaldb.h5"
#"/data020/scratch/carosko/ccha/data/ILTJ1327_16sec_1chan_testsky.h5"
f_h5 = h5parm(name_h5 ,readonly=True)
tectab = f_h5.getSoltab('sol000','tec000')
clocktab = f_h5.getSoltab('sol000','clock000')
    ####/!\ ignoring weight
freq=(f_h5.getSoltab('sol000','amplitude000').freq[0]+f_h5.getSoltab('sol000','amplitude000').freq[-1])/2#129685974.12 ###Hz taking a central freq just for the conversion of tec to sec
scal=59.*pow((150000000./freq),2)

if not(tectab.time[:]==clocktab.time[:]).all():
     ### just to make sure the time steps are the same between the different tables
     ### But no break; just information of the user from now on
    print " /!\!!!!!!!!!!!!!!!!!!!! the time steps are not the same between the tec and the clock tables"
    
nbr_tps=len(tectab.time[:])
min_tps=min(tectab.time[:])
pol=0### 1st polarisation
#times steps:
time_h5=range(nbr_tps)
#time in sec:
#for i in range(0,nbr_tps):
    #time_h5.append(tectab.time[i]-min_tps)





### A spot of house keeping
stations = ['RS106', 'RS205', 'RS208', 'RS210', 'RS305', 'RS306', 'RS307', 'RS310', 'RS406', 'RS407', 'RS409', 'RS503', 'RS508', 'RS509', 'DE601', 'DE602', 'DE603', 'DE604', 'DE605', 'FR606', 'SE607', 'UK608', 'ST001']

pdfs = []
#merged = PdfFileMerger()
count = 0
time_axis119 = [] #range(119)
time_axis89 = []

for i in range(89):
    #time_axis89.append(i * 1.337)
    time_axis89.append(i * nbr_tps/89.)
for i in range(119):
    #time_axis89.append(i * 1.337)
    time_axis119.append(i * nbr_tps/119.)    
    
    
### Make the plot
for s in stations:

    ''' CASA --------------------------------------------------------------------------------------------------- '''

    casa_data_time = []
    casa_data_pol1 = []
    casa_data_pol2 = []

    for j in range(len(x)):
        if int(x[j]) == count:
            casa_data_time.append(w[j])
            casa_data_pol1.append(y[j])
            casa_data_pol2.append(z[j])

    count = count + 1

    for i in range(len(casa_data_pol1)):
        if casa_data_pol1[i] > 350 or casa_data_pol1[i] < -350:
            casa_data_pol1[i] = casa_data_pol1[i - 1]
            
    ''' AIPS --------------------------------------------------------------------------------------------------- '''

    station_time = []
    station_pol1 = []
    station_pol2 = []
    station = s

    for i in range(len(time)):
        if frequency[i] == 129685974.12 and antenna[i] == station:
            station_time.append(time[i])
            station_pol1.append(pol1[i] * 1e9)

    for i in range(len(station_pol1)):
        if station_pol1[i] > 350 or station_pol1[i] < -200:
            station_pol1[i] = station_pol1[i - 1]
            
    ''' LOSOTO ------------------------------------------------------------------------------------------------- '''


    
    sinh5=False
    if (s+"HBA") in tectab.ant[:]:
        sinh5=True
        ant=np.where(tectab.ant[:]==s+"HBA")[0][0]### very ugly and quick solution to match the station
        

    
    ''' Plot --------------------------------------------------------------------------------------------------- '''

    plt.figure(figsize = (10, 8))
    plt.title('Polarisation 1 for ' + station)
    plt.xlabel('Time')
    plt.ylabel('Delay (ns)')

    plt.plot(time_axis89, casa_data_pol1, color = 'red', label = 'CASA')    
    plt.plot(time_axis119, station_pol1, label = 'AIPS')
    if sinh5:
       

        
        if not all(vl==0 for vl in clocktab.val[:,ant,pol]) and not all(vl==0 for vl in tectab.val[:,ant,pol]):
            plt.plot(time_h5, np.where(clocktab.val[:,ant,pol]!=0,1.e9*clocktab.val[:,ant,pol],float('NaN'))+np.where(tectab.val[:,ant,pol]!=0,scal*tectab.val[:,ant,pol],float('NaN')), label = 'Losoto')
        else: ### RQ: on pwai juste faire des jolies 0 partt)
            plt.plot(time_h5, 1.e9*clocktab.val[:,ant,pol]+scal*tectab.val[:,ant,pol], label = 'Losoto')
        
    
    plt.legend()
    plt.savefig('station-' + s + '.pdf', bbox_inches = 'tight', format = 'pdf')
    pdfs.append('station-' + s + '.pdf')
    plt.clf()


### Sort of the pdfs
output_pdf_name="delay_a.pdf"
make_pdf = "gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r150 -sOutputFile="+output_pdf_name

for s in stations:
    pdf_i='station-' + s + '.pdf'
    make_pdf=make_pdf+" "+pdf_i+" "
os.system(make_pdf)
for pdf in pdfs:
    os.remove(pdf)

transfer = "curl --upload-file ./"+output_pdf_name+" https://transfer.sh/"+output_pdf_name
os.system(transfer)
print


