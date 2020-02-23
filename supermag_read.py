#!/usr/bin/env python
#
# read supermag style csv data
#
import numpy as n
import matplotlib.pyplot as plt
import datetime
import time
import pandas
import supermag_stations as ss

class supermag_data:
    def __init__(self,fname="test_bl.csv",plot=False):
        sms = ss.supermag_stations()
        
        f=open(fname,"r")
        stats=[]
        dtimes=[]
        neds=[]
        lats=[]
        lons=[]
        for li,l in enumerate(f.readlines()):
            if li > 0:
                line=l.split(",")
                date=line[0]
                d=datetime.datetime.strptime("%sZ"%(date), "%Y-%m-%d %H:%M:%SZ")
                unixtime = time.mktime(d.timetuple())
                dtimes.append(unixtime)
                #        print(d)
                #       print(unixtime)
                stat=line[1]

                p=sms.prop(stat)
                lats.append(p["glat"])
                lons.append(p["glon"])
                stats.append(stat)
                ned=n.array([float(line[6]),float(line[7]),float(line[8])])
                neds.append(ned)
                #        print(ned)
                #       print(li)
                #      print(l)
                #stat_names=n.unique(stats)
        dtimes=n.array(dtimes)
        neds=n.array(neds)
        stats=n.array(stats)
        self.lons=n.array(lons)
        self.lats=n.array(lats)
        if plot:
            idx=n.where(stats == 'ABK')[0]
            plt.plot(dtimes[idx],neds[idx,0])
            plt.plot(dtimes[idx],neds[idx,1])
            plt.plot(dtimes[idx],neds[idx,2])
            plt.show()
            print(idx)
        self.times=dtimes
        self.ned=neds
        self.stats=stats

        
    def get_bounds(self):
        return([n.min(self.times),n.max(self.times)])

    def get_meas(self,t0,t1):
        idx=n.where((self.times >= t0)&(self.times <= t1))[0]
        return({"times":self.times[idx],
                "ned":self.ned[idx,:],
                "stat":self.stats[idx],
                "glat":self.lats[idx],
                "glon":self.lons[idx]})

if __name__ == "__main__":
    d=supermag_data(plot=False)
    b=d.get_bounds()
    print(b)
    for ti in range(24*2):
        m=d.get_meas(b[0]+ti*3600*0.5,b[0]+59+ti*3600*0.5)
        print(m)
        plt.plot(n.repeat(ti,len(m["glat"])),50*m["glat"]+m["ned"][:,0],".")

        #    plt.plot(m["ned"][:,1])
        #   plt.plot(m["ned"][:,2])
    plt.show()
