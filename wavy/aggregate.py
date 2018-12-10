from satmod import sentinel_altimeter as sa
from datetime import datetime, timedelta
#sdate=datetime(2018,8,1)
#edate=datetime(2018,9,10)
now=datetime.now()
past=datetime.now()-timedelta(days=1)
sdate=datetime(past.year,past.month,past.day)
edate=datetime(now.year,now.month,now.day)
tmpdate=sdate
timewin=0
outpath='/lustre/storeA/project/fou/om/altimeter/monthly/'
#outpath='/lustre/storeA/project/fou/om/altimeter/daily/'
while tmpdate<edate:
    tmpstart = tmpdate
    tmpend = tmpdate+timedelta(days=1)
    print ('processing: ' + str(tmpdate) + ' - ' + str(tmpend))
    sa_obj = sa(tmpstart,edate=tmpend,timewin=timewin,region="ARCMFC")
    #sa_obj = sa(tmpstart,edate=tmpend,timewin=timewin,region="Global")
    #sa_obj.dumptonc("outpath/",ncmode='auto'
    sa_obj.dumptonc(outpath,ncmode='auto',timeframe='monthly')
    tmpdate = tmpend

print ('# --- Finished --- #')
