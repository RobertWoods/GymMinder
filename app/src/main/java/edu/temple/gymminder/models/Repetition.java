package edu.temple.gymminder.models;

import com.fastdtw.timeseries.TimeSeries;

/**
 * Created by rober_000 on 4/30/2017.
 */

public class Repetition {
    int s, e, peakIndex;
    float peakValue;
    double[] accelerationStream;

    public Repetition(int s, int e, int peakIndex, float peakValue, TimeSeries t1){
        this.s = s;
        this.e = e;
        this.peakIndex = peakIndex;
        this.peakValue = peakValue;
        accelerationStream = new double[t1.size()];
        for(int i=0; i<t1.size();i++){
            accelerationStream[i] = t1.getMeasurement(i, 0);
        }
    }
}
