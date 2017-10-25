import os
import numpy as np
import matplotlib.pyplot as plt
from losoto.h5parm import h5parm

### Loading data
# cf the a for only the common stations
form="pdf"

# h5 close to aips
name_ha="/data020/scratch/lb_bw/test_losoto_delay_calibration/16sec_1SB_2comp_uvcutCS/ILTJ1327_16sec_1SB_2comp_uvcutCS_globaldb.h5"
f_ha = h5parm(name_ha ,readonly=True)
tectab_a = f_ha.getSoltab('sol000','tec000')
clocktab_a = f_ha.getSoltab('sol000','clock000')
if len(tectab_a.time[:])!=len(clocktab_a.time[:]):
    print "WARNING: The TEC and the CLOCK table don't have the same number of times points (might lead to wrong computation)"

# previous h5
name_h5="/data020/scratch/lb_bw/test_losoto_delay_calibration/noCS_circular/ILTJ1327_noCS_globaldb.h5"
f_h5 = h5parm(name_h5 ,readonly=True)
tectab_5 = f_h5.getSoltab('sol000','tec000')
clocktab_5 = f_h5.getSoltab('sol000','clock000')
if len(tectab_5.time[:])!=len(clocktab_5.time[:]):
    print "WARNING: The TEC and the CLOCK table don't have the same number of times points (might lead to wrong computation)"


#~ freq=129685974.12 ###MHz taking a central freq just for the conversion of tec to sec
#~ scal=59.*pow((150000000./freq),2)


n_ant=-1
for ant in tectab_a.ant[:]:
    n_pol=-1
    n_ant=n_ant+1
    print tectab_a.ant[n_ant]
    sinh5=False
    if ant in tectab_5.ant[:]:
        sinh5=True
        n_ant_5=np.where(tectab_5.ant[:]==ant)[0][0]### very ugly and quick solution to match the station
        
    
    
    print ant
    fig = plt.figure()
    fig.suptitle(ant)
    
    for pol in tectab_a.pol[:]:
        n_pol=n_pol+1
       
        if n_pol==0:
            ax_tecxx = fig.add_subplot(221)
            plt.ylabel('TEC ($10^{16} \cdot \mathrm{m}^{-2}$)')
            plt.setp(ax_tecxx.get_xticklabels(), visible=False)
            
        else:
            ax_tecyy = fig.add_subplot(222,sharey=ax_tecxx)
            plt.setp(ax_tecyy.get_yticklabels(), visible=False)
            plt.setp(ax_tecyy.get_xticklabels(), visible=False)
            
        ### Assuming all 0 value are Nan...
        if not all(vl==0 for vl in tectab_a.val[:,n_ant,n_pol]): ### just to not erase completely the reference antenna
            plt.plot(clocktab_a.time[:]-clocktab_a.time[0],np.where(tectab_a.val[:,n_ant,n_pol]!=0,tectab_a.val[:,n_ant,n_pol],float('Nan')), label="H5 previous")
        else:
            plt.plot(clocktab_a.time[:]-clocktab_a.time[0],tectab_a.val[:,n_ant,n_pol])### assuming a file has always the same colors , label="H5 previous")
        
        if sinh5:
            if not all(vl==0 for vl in tectab_5.val[:,n_ant_5,n_pol]):
                plt.plot(clocktab_a.time[:]-clocktab_a.time[0],np.where(tectab_5.val[:,n_ant_5,n_pol]!=0,tectab_5.val[:,n_ant_5,n_pol],float('NaN')), label="H5 close to AIPS")
            else:
                plt.plot(clocktab_a.time[:]-clocktab_a.time[0],tectab_5.val[:,n_ant_5,n_pol])###, label="H5 close to AIPS")
        plt.title(pol)
        
        
        
        if n_pol==0:
            ax_clockxx = fig.add_subplot(223,sharex=ax_tecxx)
            plt.ylabel('Clock (ns)')
        else:
            ax_clockyy = fig.add_subplot(224,sharex=ax_tecxx,sharey=ax_clockxx) acquis
            plt.setp(ax_clockyy.get_yticklabels(), visible=False)
        
        
        if not all(vl==0 for vl in clocktab_a.val[:,n_ant,n_pol]): ### just to not erase completely the reference antenna
            if n_pol==1:
                plt.plot(clocktab_a.time[:]-clocktab_a.time[0],np.where(clocktab_a.val[:,n_ant,n_pol]!=0,clocktab_a.val[:,n_ant,n_pol]*1.e9,float('Nan')), label="H5 previous")
            else:
                plt.plot(clocktab_a.time[:]-clocktab_a.time[0],np.where(clocktab_a.val[:,n_ant,n_pol]!=0,clocktab_a.val[:,n_ant,n_pol]*1.e9,float('Nan')))###, label="H5 previous")
        else:
            if n_pol==1:
                plt.plot(clocktab_a.time[:]-clocktab_a.time[0],clocktab_a.val[:,n_ant,n_pol]*1.e9, label="H5 previous")
            else:
                plt.plot(clocktab_a.time[:]-clocktab_a.time[0],clocktab_a.val[:,n_ant,n_pol]*1.e9)###, label="H5 previous")
                
        if sinh5:
            if not all(vl==0 for vl in clocktab_5.val[:,n_ant_5,n_pol]):
                if n_pol==1:
                    plt.plot(clocktab_a.time[:]-clocktab_a.time[0],np.where(clocktab_5.val[:,n_ant_5,n_pol]!=0,clocktab_5.val[:,n_ant_5,n_pol]*1.e9,float('NaN')), label="H5 close to AIPS")
                else:
                    plt.plot(clocktab_a.time[:]-clocktab_a.time[0],np.where(clocktab_5.val[:,n_ant_5,n_pol]!=0,clocktab_5.val[:,n_ant_5,n_pol]*1.e9,float('NaN')))###, label="H5 close to AIPS")
            else:
                if n_pol==1:
                    plt.plot(clocktab_a.time[:]-clocktab_a.time[0],clocktab_5.val[:,n_ant_5,n_pol]*1.e9, label="H5 close to AIPS")
                else:
                    plt.plot(clocktab_a.time[:]-clocktab_a.time[0],clocktab_5.val[:,n_ant_5,n_pol]*1.e9)###, label="H5 close to AIPS")

        plt.xlabel('Time (s)')
        plt.legend()
    
    
       
    plt.savefig('station-' + ant +'.'+form, bbox_inches = 'tight', format = form)
    plt.clf()


#### Sort of the pdfs
output_pdf_name="clocktec_h5s."+form
make_pdf = "gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r150 -sOutputFile="+output_pdf_name
pdfs = []
for ant in tectab_a.ant[:]:
    for pol in tectab_a.pol[:]:
        pdf_i='station-' + ant +'.'+form
        #pdfs.append(pdf_i)
        make_pdf=make_pdf+" "+pdf_i+" "
os.system(make_pdf)
for pdf in pdfs:
    os.remove(pdf)

transfer = "curl --upload-file ./"+output_pdf_name+" https://transfer.sh/"+output_pdf_name
os.system(transfer)
print

        
