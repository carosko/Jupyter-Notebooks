#!/usr/bin/env python




# import required modules ---------------------------------------------------------------------

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from losoto.h5parm import h5parm


def input2strlist_nomapfile(invar):    
    """ from (prefactor) bin/download_IONEX.py
    give the list of MSs/files from the list provided as a string
    """
   
    str_list = None
    if type(invar) is str:
        if invar.startswith('[') and invar.endswith(']'):
            str_list = [f.strip(' \'\"') for f in invar.strip('[]').split(',')]
        else:
            str_list = [invar.strip(' \'\"')]
    elif type(invar) is list:
        str_list = [str(f).strip(' \'\"') for f in invar]
    else:
        raise TypeError('input2strlist: Type '+str(type(invar))+' unknown!')
    return str_list


def clockh5_plot(h5_file,res_dir):
    ''' clock versus time for XX
        AND
        clock difference between XX and YY versus time
    '''
    
    # load the data -------------------------------------------------------------------------------

    #name = "/data020/scratch/lb_bw/test_losoto_delay_calibration/noCS_circular/ILTJ1327_noCS_globaldb.h5" #'/home/sean/test.h5'
    #name="/data020/scratch/lb_bw/test_losoto_delay_calibration/16sec_1SB_2comp_uvcutCS/ILTJ1327_16sec_1SB_2comp_uvcutCS_globaldb.h5"
    f = h5parm(h5_file, readonly = True)
    clock = f.getSoltab('sol000', 'clock000')


    ''' clock versus time for XX ============================================================== '''

    # set up the plot -----------------------------------------------------------------------------

    p = 0 # 0 = XX; 1 = YY

    plt.figure(figsize = (24, 12))
    normalised = [(time - min(clock.time[:]))/(max(clock.time[:]) - min(clock.time[:])) for time in clock.time[:]]


    # plot the clock values for each station ------------------------------------------------------

    for a in (range(len(clock.ant[:]))):
      
        ### About the title and the "crop" of the names of stations: maybe I am not using the same data file, but it doesn't seem to be needed for me :/... (would erase the letters at the beginning of the name actually). Same for the polarization
       
        ### (not ignoring the "jump" to zero
        #plt.plot(normalised, 1.e9*clock.val[:,a,p], lw = 0.5, label = clock.ant[a])
       
        ### assuming that ALL zeros are previous NaN
        if not all(vl==0 for vl in clock.val[:,a,p]):### just to not flagged completely the reference antenna
            plt.plot(normalised, np.where(clock.val[:,a,p]!=0,1.e9*clock.val[:,a,p],float('NaN')), lw = 0.5, label = clock.ant[a])
        else:
            plt.plot(normalised, 1.e9*clock.val[:,a,p], lw = 0.5, label = clock.ant[a])
        

    # tidy up the plot ----------------------------------------------------------------------------

    mpl.rcParams.update({'font.size': 10})
    plt.xlim(min(normalised), max(normalised))
    ###plt.title('Clock versus time for the ' + str(clock.pol[p])[2:-1] + ' polarisation for all stations', fontsize = 10)
    plt.title('Clock versus time for the ' + clock.pol[p] + ' polarisation for all stations', fontsize = 10)
    plt.xlabel('Time (normalised)', fontsize = 10)
    plt.ylabel(r'Clock (ns)', fontsize = 10)
    plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.05), ncol = 17, frameon = False, fontsize = 7)
    plt.savefig(res_dir+'/clock_xx.png', bbox_inches = 'tight')
    plt.close('all')



    ''' clock difference between XX and YY versus time ======================================== '''

    # set up the plot -----------------------------------------------------------------------------

    plt.figure(figsize = (24, 12)) ###[needed to be repeated?]
    #~ normalised = [(time - min(clock.time[:]))/(max(clock.time[:]) - min(clock.time[:])) for time in clock.time[:]]

    # plot the clock difference for each station --------------------------------------------------

    for a in (range(len(clock.ant[:]))):
       
        ### np.where or the plot might be able to complain that xx and yy have not the same lenght -or at least manage with that
        ### but just to be sure, let's keep it:
        if len(clock.val[:,a,0])!=len(clock.val[:,a,1]):
            print('Cannot find difference as lists of XX and YY clock values have different lengths. Exiting...')
            sys.exit()
        
      
        
        ### with no 0/NaN flagging
        ###plt.plot(normalised, 1.e9*(clock.val[:,a,0]-clock.val[:,a,1]), lw = 0.5, label = clock.ant[a])
      
        ### with the 0/Nan flagging
        if (not all(vl==0 for vl in clock.val[:,a,0])) and (not all(vl==0 for vl in clock.val[:,a,1])) :
            ### (NB: We want to flag the 0/NaN before making the difference
            ###     as we are expecting the difference between the polarization to be very small or most of the time equal to 0)
            plt.plot(normalised, np.where(clock.val[:,a,0]!=0,1.e9*clock.val[:,a,0],float('NaN'))-np.where(clock.val[:,a,1]!=0,1.e9*clock.val[:,a,1],float('NaN')), lw = 0.5, label = clock.ant[a])
        else:
            plt.plot(normalised, 1.e9*(clock.val[:,a,0]-clock.val[:,a,1]), lw = 0.5, label = clock.ant[a])

    # tidy up the plot ----------------------------------------------------------------------------

    mpl.rcParams.update({'font.size': 10})
    plt.xlim(min(normalised), max(normalised))
    #plt.title('Polarisation ' + str(clock.pol[0])[2:-1] + ' and ' + str(clock.pol[1])[2:-1] + ' clock difference versus time for all stations', fontsize = 10)
    plt.title('Polarisation ' + clock.pol[0] + ' and ' + clock.pol[1] + ' clock difference versus time for all stations', fontsize = 10)
    plt.xlabel('Time (normalised)', fontsize = 10)
    #plt.ylabel(str(clock.pol[0])[2:-1] + ' and ' + str(clock.pol[1])[2:-1] + ' clock difference (ns)', fontsize = 10)
    plt.ylabel(clock.pol[0] + ' and ' + clock.pol[1] + ' clock difference (ns)', fontsize = 10)
    plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.05), ncol = 17, frameon = False, fontsize = 7)
    plt.savefig(res_dir+'/clock_difference_xx_yy.png', bbox_inches = 'tight')
    plt.close('all')



def tech5_plot(h5_file,res_dir):
    ''' TEC versus time for XX ================================================================ 
        AND
        TEC difference between XX and YY versus time
    '''


    # load the data -------------------------------------------------------------------------------
    f = h5parm(h5_file, readonly = True)
    tec = f.getSoltab('sol000', 'tec000')
    
    
    ''' TEC versus time for XX ================================================================ '''
    
    # set up the plot -----------------------------------------------------------------------------

    p = 0 # 0 = XX; 1 = YY
    plt.figure(figsize = (24, 12))
    normalised = [(time - min(tec.time[:]))/(max(tec.time[:]) - min(tec.time[:])) for time in tec.time[:]]


    # plot the TEC values for each station --------------------------------------------------------

    for a in (range(len(tec.ant[:]))):
        ### Without 0/NaN flagging
        ###plt.plot(normalised, tec.val[:,a,p], lw = 0.5, label = tec.ant[a])
        
        ### With flagging of ALL 0/NaN
        if not all(vl==0 for vl in tec.val[:,a,p]):### just to not flagged completely the reference antenna
            plt.plot(normalised, np.where(tec.val[:,a,p]!=0,tec.val[:,a,p],float('NaN')), lw = 0.5, label = tec.ant[a])
        else:
            plt.plot(normalised, tec.val[:,a,p], lw = 0.5, label = tec.ant[a])
            

    # tidy up the plot ----------------------------------------------------------------------------

    mpl.rcParams.update({'font.size': 10})
    plt.xlim(min(normalised), max(normalised))
    #plt.title('TEC versus time for the ' + str(tec.pol[p])[2:-1] + ' polarisation for all stations', fontsize = 10)
    plt.title('TEC versus time for the ' + tec.pol[p] + ' polarisation for all stations', fontsize = 10)
    plt.xlabel('Time (normalised)', fontsize = 10)
    plt.ylabel(r'TEC ($10^{16} \cdot \mathrm{m}^{-2}$)', fontsize = 10)
    plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.05), ncol = 17, frameon = False, fontsize = 7)
    plt.savefig(res_dir+'/tec_xx.png', bbox_inches = 'tight')
    plt.close('all')




    ''' TEC difference between XX and YY versus time ========================================== '''

    # set up the plot -----------------------------------------------------------------------------

    plt.figure(figsize = (24, 12))
    #~ normalised = [(time - min(tec.time[:]))/(max(tec.time[:]) - min(tec.time[:])) for time in tec.time[:]]


    # plot the TEC difference for each station ----------------------------------------------------

    for a in (range(len(tec.ant[:]))):
        
        if len(tec.val[:,a,0])!=len(tec.val[:,a,1]):
            print('Cannot find difference as lists of XX and YY clock values have different lengths. Exiting...')
            sys.exit()
        

        ### Without 0/NaN flagging
        ###plt.plot(normalised, tec.val[:,a,0]-tec.val[:,a,1], lw = 0.5, label = tec.ant[a])
        
        ### with the 0/Nan flagging
        if (not all(vl==0 for vl in tec.val[:,a,0])) and (not all(vl==0 for vl in tec.val[:,a,1])) :
            plt.plot(normalised, np.where(tec.val[:,a,0]!=0,tec.val[:,a,0],float('NaN'))-np.where(tec.val[:,a,1]!=0,tec.val[:,a,1],float('NaN')), lw = 0.5, label = tec.ant[a])
        else:
            plt.plot(normalised, tec.val[:,a,0]-tec.val[:,a,1], lw = 0.5, label = tec.ant[a])

      
        
    # tidy up the plot ----------------------------------------------------------------------------

    mpl.rcParams.update({'font.size': 10})
    plt.xlim(min(normalised), max(normalised))
    #plt.title('Polarisation ' + str(tec.pol[0])[2:-1] + ' and ' + str(tec.pol[1])[2:-1] + ' TEC difference versus time for all stations', fontsize = 10)
    plt.title('Polarisation ' + tec.pol[0] + ' and ' + tec.pol[1] + ' TEC difference versus time for all stations', fontsize = 10)
    plt.xlabel('Time (normalised)', fontsize = 10)
    #plt.ylabel(str(tec.pol[0])[2:-1] + ' and ' + str(tec.pol[1])[2:-1] + r' TEC difference ($10^{16} \cdot \mathrm{m}^{-2}$)', fontsize = 10)
    plt.ylabel(tec.pol[0] + ' and ' + tec.pol[1] + r' TEC difference ($10^{16} \cdot \mathrm{m}^{-2}$)', fontsize = 10)
    plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.05), ncol = 17, frameon = False, fontsize = 7)
    plt.savefig(res_dir+'/tec_difference_xx_yy.png', bbox_inches = 'tight')
    plt.close('all')


def amph5_xxyy_plot(h5_file,res_dir):
    ''' amplitude versus frequency for XX and YY ============================================== '''
    

    # load the data -------------------------------------------------------------------------------
    f = h5parm(h5_file, readonly = True)
    amplitude = f.getSoltab('sol000', 'amplitude000')


    # set up the plot -----------------------------------------------------------------------------

    plt.figure(figsize = (24, 12))
    normalised = [(frequency - min(amplitude.freq[:]))/(max(amplitude.freq[:]) - min(amplitude.freq[:])) for frequency in amplitude.freq[:]]

    # plot the amplitude for each station ---------------------------------------------------------

    # note: I think the plot is correct but all of the core stations have one massive peak that
    #       is drowing out all of the other points on the plot so I added a step that clips all
    #       values above a given threshold

    threshold =1000. #1.e6 #3

    text = open(res_dir+'/clippings.txt', 'w')
    text.write('Clipping threshold: ' + str(threshold) + '\n')

    # note: here I'm looping over all of the stations and for each polarisation I end up with a
    #       list which has 843 elements - one for each time step - and each of those 843 time
    #       steps is a list itself of 122 elements - one value for each subband

    for a in range(len(amplitude.ant[:])):

    # note: the xx list has the 843 elements and the average_xx list has 122 - one for each
    #       subband which was found by averaging over the 843 time steps for each frequency
    #       value, and this is done for the X and Y polarisations
        for p in [0,1]:
           
            ### NB: we keep "average" (and not median) for now to not get rid artificialy of the outliers
                ### [TODO: double check the axis for average: seems to be good though] 
                ### (amplitude.val[polarisation, direction, station, frequency, time])
            #~ plt.plot(normalised, np.average(amplitude.val[p,0,a,:,:],axis=-1), lw = 0.5, label = amplitude.ant[a] + ' ' + amplitude.pol[p])
           
            ### (just to have the result of a median by hand)
            #~ plt.plot(normalised, np.median(amplitude.val[p,0,a,:,:],axis=-1), lw = 0.5, label = amplitude.ant[a] + ' ' + amplitude.pol[p])
           
            ### let's try to flag the pics:
            ### (the outliers will be put at the threshold value - and not to the previous value
            ###     (I don't know how to do that in one line :$ ))
            ### NB: if there are NaNs put to 0, they will lower the "true" average value
            plt.plot(normalised, np.where(np.average(amplitude.val[p,0,a,:,:],axis=-1)>threshold,threshold,np.average(amplitude.val[p,0,a,:,:],axis=-1)), lw = 0.5, label = amplitude.ant[a] + ' ' + amplitude.pol[p]) 
            
            ### and finally a small test just to check the influence of 0/NaN
            ###     (which can show that 0/NaN have naturally no influence "for now" when the pic is present 
            #~ plt.plot(normalised, np.where(np.average(np.where(amplitude.val[p,0,a,:,:]!=0,amplitude.val[p,0,a,:,:],float('NaN')),axis=-1)>threshold,threshold,np.average(np.where(amplitude.val[p,0,a,:,:]!=0,amplitude.val[p,0,a,:,:],float('NaN')),axis=-1)), lw = 0.5, label = amplitude.ant[a] + ' ' + amplitude.pol[p])
            
            ### (a way to filter the clipped values
            ### The non flagged data are put to NaN
            ### So the length of the clipped array should be the same than for the amplitude 
            ### )
            clipped_pol=[]
            clipped_pol=np.where(np.average(amplitude.val[p,0,a,:,:],axis=-1)<threshold,float('NaN'),np.average(amplitude.val[p,0,a,:,:],axis=-1))
            #~ text.write('Station: ' + str(amplitude.ant[a])[2:-1] + '; Polarisation: ' + str(amplitude.pol[p])[2:-1] + '; Clipped: ' + str(np.count_nonzero(~np.isnan(clipped_pol))) + '\n')    
            text.write('Station: ' + amplitude.ant[a] + '; Polarisation: ' + amplitude.pol[p] + '; Clipped: ' + str(np.count_nonzero(~np.isnan(clipped_pol))) + '\n')    

    #~ # note: this is where the clipping is actually done, where the values are saved to a new
    #~ #       list if they meet the criteria, and if they don't I add the most recent value that
    #~ #       did meet the criteria so the plot looks continuous (except where there is no recent
    #~ #       value and then 0 is added)
            
        #~ count_x, count_y, clipped_xx, clipped_yy = [], [], [], []

        #~ ### To "glu" the clipping with the modifications done above
        #~ p=0
        #~ average_xx=np.average(amplitude.val[p,0,a,:,:],axis=-1)
        #~ p=1
        #~ average_yy=np.average(amplitude.val[p,0,a,:,:],axis=-1)
        #~ ####
        #~ if len(average_xx) == len(average_yy):
            #~ for clip in range(len(average_xx)):
                #~ if average_xx[clip] < threshold:
                    #~ clipped_xx.append(np.average(average_xx[clip]))

                #~ else:
                    #~ if len(clipped_xx) == 0:
                        #~ clipped_xx.append(0)
                    
                    #~ else:
                        #~ clipped_xx.append(clipped_xx[clip - 1])

                    #~ count_x.append(0)

                #~ if average_yy[clip] < threshold:
                    #~ clipped_yy.append(np.average(average_yy[clip]))

                #~ else:
                    #~ if len(clipped_yy) == 0:
                        #~ clipped_yy.append(0)
                    
                    #~ else:
                        #~ clipped_yy.append(clipped_yy[clip - 1])

                    #~ count_y.append(0)

        #~ text.write('Station: ' + str(amplitude.ant[a])[2:-1] + '; Polarisation: ' + str(amplitude.pol[0])[2:-1] + '; Clipped: ' + str(len(count_x)) + '\n')
        #~ text.write('Station: ' + str(amplitude.ant[a])[2:-1] + '; Polarisation: ' + str(amplitude.pol[1])[2:-1] + '; Clipped: ' + str(len(count_y)) + '\n')

        #~ plt.plot(normalised, clipped_xx, lw = 0.5, label = str(amplitude.ant[a])[2:-1] + ' ' + str(amplitude.pol[0])[2:-1]) # one series per station per polarisation
        #~ plt.plot(normalised, clipped_yy, lw = 0.5, label = str(amplitude.ant[a])[2:-1] + ' ' + str(amplitude.pol[1])[2:-1])

    text.close()
        
    # tidy up/ the plot ----------------------------------------------------------------------------

    mpl.rcParams.update({'font.size': 10})
    plt.xlim(min(normalised), max(normalised))
    #plt.title('Amplitude versus frequency for ' + str(amplitude.pol[0])[2:-1] + ' and ' + str(amplitude.pol[1])[2:-1] + ' polarisations', fontsize = 10)
    plt.title('Amplitude versus frequency for ' + amplitude.pol[0] + ' and ' + amplitude.pol[1] + ' polarisations', fontsize = 10)
    plt.xlabel('Frequency (normalised)', fontsize = 10)
    #plt.ylabel(str(amplitude.pol[0])[2:-1] + ' and ' + str(amplitude.pol[1])[2:-1] + ' amplitude', fontsize = 10)
    plt.ylabel(amplitude.pol[0]+ ' and ' + amplitude.pol[1] + ' amplitude', fontsize = 10)
    plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.05), ncol = 17, frameon = False, fontsize = 7)
    plt.savefig(res_dir+'/amplitude.png', bbox_inches = 'tight')
    plt.close('all')


def main(h5_list, output_dir):
    
    ### extraction of one h5 for clock and tec plots
    h5_0=input2strlist_nomapfile(h5_list)[0]
    clockh5_plot(h5_0,output_dir)
    tech5_plot(h5_0,output_dir)
    
    ### TODO for amplitudes: we might want to plot all SBs at once (if available)
    amph5_xxyy_plot(h5_0,output_dir)


################
#main("/data020/scratch/lb_bw/test_losoto_delay_calibration/16sec_1SB_2comp_uvcutCS/ILTJ1327_16sec_1SB_2comp_uvcutCS_globaldb.h5","./res")

########################################################################
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Plot the solution for the h5 file provided\n\t(amp VS freq, clock and tec (and their differences between polarisations) VS time, ')
    
    parser.add_argument('h5', type=str, nargs='+',
                        help='One (or more) h5(s) for which the plots are done')
    parser.add_argument('--ResDir', type=str, 
                        help='Path to where to store the outputs (plots and count of the clipped points)')
    
   
    args = parser.parse_args()
    res_dir="./"
    if args.ResDir:
        res_dir=args.ResDir
    
    main(args.h5,res_dir)

