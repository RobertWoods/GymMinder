package edu.temple.gymminder;

import com.fastdtw.timeseries.TimeSeriesBase;

import org.junit.rules.ExternalResource;

import java.util.ArrayList;

import edu.temple.gymminder.processing.DataUtils;

/**
 * Created by rober_000 on 4/12/2017.
 */

public class DataUtilsResources extends ExternalResource {

    public final double PERIOD = 1/10000L;
    public final int S2MS_CONVERSION = 1000000;
    ArrayList<ArrayList<Float>> data = new ArrayList<>();
    ArrayList<ArrayList<Float>> processed = new ArrayList<>();
    ArrayList<Long> timestamps;

    protected void before(){
        for (int i = 0; i < 3; i++) {
            data.add(new ArrayList<Float>());
            processed.add(new ArrayList<Float>());
        }
        timestamps = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            timestamps.add((long) PERIOD * i * S2MS_CONVERSION);
        }
        DataUtils.init(data, timestamps, processed);
    }

    public Class getAccessibleDataUtilsClass(String className) throws NoSuchFieldException {
        Class[] classes = DataUtils.class.getDeclaredClasses();
        Class peakClass = null;
        for (Class c : classes) {
            if (c.getName().equals("edu.temple.gymminder.processing.DataUtils$" + className)) {
                peakClass = c;
                break;
            }
        }
        switch(className){
            case "Peak":
                peakClass.getDeclaredField("amplitude").setAccessible(true);
                peakClass.getDeclaredField("index").setAccessible(true);
                break;
            case "DetectedBounds":
                //Don't need any fields for now o3o
                break;
        }
        return peakClass;
    }

    void setupPeakTimeSeriesAndAxis(){
        TimeSeriesBase.Builder builder = new TimeSeriesBase.Builder();
        for(int i=0;i<50;i++){
            builder = builder.add(0, i);
        }
        DataUtils.repTimeSeries = builder.build();
        DataUtils.repPeak = new DataUtils.Peak(25,0);
        DataUtils.majorAxisIndex = 0;
        DataUtils.setListener(null);
        DataUtils.init(data, timestamps, processed);
    }

}
