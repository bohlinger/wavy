#!/usr/bin/env python3
import os
import sys
from datetime import datetime, timedelta

now = datetime(2019,12,22,6)
init_fname = now.strftime("%Y%m%dT%H")
outpath = "/lustre/storeA/project/fou/om/LC-WVF/"

lt_lst = ['00','03','06','09','12','15','18','21','24','27','30','33','36','39','42','45','48','51','54','57','60']

for lt in lt_lst:
    print('leadtime: ', lt)
    fc_str = (now + timedelta(hours=int(lt))).strftime("%Y-%m-%dT%H:00:00")
    fc_fname = (now + timedelta(hours=int(lt))).strftime("%Y%m%dT%H" )

    ustr = ("fimex "
            + "--input.file=/lustre/storeB/immutable/archive/projects/metproduction/DNMI_WAVE/" 
            + now.strftime("%Y") 
            + "/" 
            + now.strftime("%m") 
            + "/" 
            + now.strftime("%d") 
            + "/W4km_force_" 
            + init_fname + "Z.nc"
            + " --input.type=nc"
            + " --interpolate.method=bilinear"
            + ' --interpolate.projString="+proj=latlong +ellps=sphere +a=6371000 +e=0"'
            + " --interpolate.xAxisValues=-32.0,-31.95,...,95.00"
            + " --interpolate.yAxisValues=40.0,40.03,...,85.00"
            + " --interpolate.xAxisUnit=degree"
            + " --interpolate.yAxisUnit=degree"
            + " --process.rotateVector.all"
            + " --extract.selectVariables=Uwind"
            + " --extract.selectVariables=forecast_reference_time"
            + " --extract.reduceTime.start=" + fc_str
            + " --extract.reduceTime.end=" + fc_str
            + " --output.config=/home/patrikb/wavy/LC-WVF/cdmGribWriterConfig.xml"
            + " --output.file=" + outpath
            + "wave_emni_" + fc_fname + "_test_fc_10u.grib2"
            + " --output.type=grib2 --input.config=/home/patrikb/wavy/LC-WVF/wind_height.ncml")
    tmp=os.system(ustr)

    vstr = ("fimex "
            + "--input.file=/lustre/storeB/immutable/archive/projects/metproduction/DNMI_WAVE/"
            + now.strftime("%Y")
            + "/"
            + now.strftime("%m")
            + "/"
            + now.strftime("%d")
            + "/W4km_force_"
            + init_fname + "Z.nc"
            + " --input.type=nc"
            + " --interpolate.method=bilinear"
            + ' --interpolate.projString="+proj=latlong +ellps=sphere +a=6371000 +e=0"'
            + " --interpolate.xAxisValues=-32.0,-31.95,...,95.00"
            + " --interpolate.yAxisValues=40.0,40.03,...,85.00"
            + " --interpolate.xAxisUnit=degree"
            + " --interpolate.yAxisUnit=degree"
            + " --process.rotateVector.all"
            + " --extract.selectVariables=Vwind"
            + " --extract.selectVariables=forecast_reference_time"
            + " --extract.reduceTime.start=" + fc_str
            + " --extract.reduceTime.end=" + fc_str
            + " --output.config=/home/patrikb/wavy/LC-WVF/cdmGribWriterConfig.xml"
            + " --output.file=" + outpath
            + "wave_emni_" + fc_fname + "_test_fc_10v.grib2"
            + " --output.type=grib2 --input.config=/home/patrikb/wavy/LC-WVF/wind_height.ncml")
    tmp=os.system(vstr)

    mp2str = ("fimex "
            + " --input.file=/lustre/storeB/immutable/archive/projects/metproduction/DNMI_WAVE/"
            + now.strftime("%Y")
            + "/"
            + now.strftime("%m")
            + "/"
            + now.strftime("%d")
            + "/MyWave_wam4_WAVE_"
            + init_fname + "Z.nc"
            + " --input.type=nc"
            + " --interpolate.method=bilinear"
            + ' --interpolate.projString="+proj=latlong +ellps=sphere +a=6371000 +e=0"'
            + " --interpolate.xAxisValues=-32.0,-31.95,...,95.00"
            + " --interpolate.yAxisValues=40.0,40.03,...,85.00"
            + " --interpolate.xAxisUnit=degree"
            + " --interpolate.yAxisUnit=degree"
            + " --process.rotateVector.all"
            + " --extract.selectVariables=tm2"
            + " --extract.selectVariables=forecast_reference_time"
            + " --extract.reduceTime.start=" + fc_str
            + " --extract.reduceTime.end=" + fc_str
            + " --output.config=/home/patrikb/wavy/LC-WVF/cdmGribWriterConfig.xml"
            + " --output.file=" + outpath
            + "wave_emni_" + fc_fname + "_test_fc_mp2.grib2" 
            + " --output.type=grib2;")
    tmp=os.system(mp2str)

    mwdstr = ("fimex "
            + " --input.file=/lustre/storeB/immutable/archive/projects/metproduction/DNMI_WAVE/"
            + now.strftime("%Y")
            + "/"
            + now.strftime("%m")
            + "/"
            + now.strftime("%d")
            + "/MyWave_wam4_WAVE_"
            + init_fname + "Z.nc"
            + " --input.type=nc"
            + " --interpolate.method=bilinear"
            + ' --interpolate.projString="+proj=latlong +ellps=sphere +a=6371000 +e=0"'
            + " --interpolate.xAxisValues=-32.0,-31.95,...,95.00"
            + " --interpolate.yAxisValues=40.0,40.03,...,85.00"
            + " --interpolate.xAxisUnit=degree"
            + " --interpolate.yAxisUnit=degree"
            + " --process.rotateVector.all"
            + " --extract.selectVariables=thq"
            + " --extract.selectVariables=forecast_reference_time"
            + " --extract.reduceTime.start=" + fc_str
            + " --extract.reduceTime.end=" + fc_str
            + " --output.config=/home/patrikb/wavy/LC-WVF/cdmGribWriterConfig.xml"
            + " --output.file=" + outpath
            + "wave_emni_" + fc_fname + "_test_fc_mwd.grib2"
            + " --output.type=grib2;")
    tmp=os.system(mwdstr)

    tpstr = ("fimex "
            + " --input.file=/lustre/storeB/immutable/archive/projects/metproduction/DNMI_WAVE/"
            + now.strftime("%Y")
            + "/"
            + now.strftime("%m")
            + "/"
            + now.strftime("%d")
            + "/MyWave_wam4_WAVE_"
            + init_fname + "Z.nc"
            + " --input.type=nc"
            + " --interpolate.method=bilinear"
            + ' --interpolate.projString="+proj=latlong +ellps=sphere +a=6371000 +e=0"'
            + " --interpolate.xAxisValues=-32.0,-31.95,...,95.00"
            + " --interpolate.yAxisValues=40.0,40.03,...,85.00"
            + " --interpolate.xAxisUnit=degree"
            + " --interpolate.yAxisUnit=degree"
            + " --process.rotateVector.all"
            + " --extract.selectVariables=tp"
            + " --extract.selectVariables=forecast_reference_time"
            + " --extract.reduceTime.start=" + fc_str
            + " --extract.reduceTime.end=" + fc_str
            + " --output.config=/home/patrikb/wavy/LC-WVF/cdmGribWriterConfig.xml"
            + " --output.file=" + outpath
            + "wave_emni_" + fc_fname + "_test_fc_tp.grib2"
            + " --output.type=grib2;")
    tmp=os.system(tpstr)

    swhstr = ("fimex "
            + " --input.file=/lustre/storeB/immutable/archive/projects/metproduction/DNMI_WAVE/"
            + now.strftime("%Y")
            + "/"
            + now.strftime("%m")
            + "/"
            + now.strftime("%d")
            + "/MyWave_wam4_WAVE_"
            + init_fname + "Z.nc"
            + " --input.type=nc"
            + " --interpolate.method=bilinear"
            + ' --interpolate.projString="+proj=latlong +ellps=sphere +a=6371000 +e=0"'
            + " --interpolate.xAxisValues=-32.0,-31.95,...,95.00"
            + " --interpolate.yAxisValues=40.0,40.03,...,85.00"
            + " --interpolate.xAxisUnit=degree"
            + " --interpolate.yAxisUnit=degree"
            + " --process.rotateVector.all"
            + " --extract.selectVariables=hs"
            + " --extract.selectVariables=forecast_reference_time"
            + " --extract.reduceTime.start=" + fc_str
            + " --extract.reduceTime.end=" + fc_str
            + " --output.config=/home/patrikb/wavy/LC-WVF/cdmGribWriterConfig.xml"
            + " --output.file=" + outpath
            + "wave_emni_" + fc_fname + "_test_fc_hs.grib2"
            + " --output.type=grib2;")
    tmp=os.system(swhstr)

cmd = ("cat " + outpath + "*test*.grib2 > " 
    + outpath + "wave_enmi_" + now.strftime("%Y%m%d%H") + "_prod_fc.grib2")
tmp=os.system(cmd)
cmd = ("rm " + outpath + "*test*")
tmp=os.system(cmd)
