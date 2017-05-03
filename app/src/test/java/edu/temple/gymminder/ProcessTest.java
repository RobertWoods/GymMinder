package edu.temple.gymminder;

import android.util.SparseArray;

import com.fastdtw.timeseries.TimeSeriesBase;

import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

/**
 * Created by rober_000 on 4/12/2017.
 */

public class ProcessTest {
    @Rule
    public final DataUtilsResources res = new DataUtilsResources();

    @Before
    public void setup(){
        res.timestamps = new ArrayList<>();
        DataUtils.init(res.data, res.timestamps, res.processed);
    }

    @Test
    public void processDoesNotCrash(){
        setupPeakTimeSeriesAndAxis();
        Random random = new Random();
        for(int i=0; i<100; i++){
            float[] values = {random.nextFloat(), random.nextFloat(), random.nextFloat()};
            DataUtils.process(values, DataUtils.POLLING_RATE*i);
        }
    }

    @Test
    public void processCorrectlyAddsValuesWithTimestampsGreaterThanPeriod() throws InterruptedException {
        processDoesNotCrash();
        assertEquals(100, res.timestamps.size());
        assertEquals(100, res.data.get(0).size());
        assertEquals(100, res.processed.get(0).size());
    }

    @Test
    public void processCorrectlyAddsValuesWithTimestampsLessThanPeriod() throws InterruptedException {
        setupPeakTimeSeriesAndAxis();
        Random random = new Random();
        for(int i=0; i<100; i++){
            float[] values = {random.nextFloat(), random.nextFloat(), random.nextFloat()};
            DataUtils.process(values, DataUtils.POLLING_RATE/2*i);
        }
        //Should be 1+floor(99/2)
        //Or 1+floor(100/2), it honestly doesn't matter
        assertEquals(50, res.timestamps.size(), 1);
        assertEquals(50, res.data.get(0).size(), 1);
        assertEquals(50, res.processed.get(0).size(), 1);
    }

    @Test
    public void processCorrectlyAddsPeak(){
        setupPeakTimeSeriesAndAxis();
        int i;
        for(i=0; i<50; i++){
            float[] values = {0, 0, 0};
            DataUtils.process(values, DataUtils.POLLING_RATE*i);
        }
        DataUtils.process(new float[] {100,100,100}, DataUtils.POLLING_RATE*i++);
        int j = i;
        for(; i<j+5; i++){
            float[] values = {0, 0, 0};
            DataUtils.process(values, DataUtils.POLLING_RATE*(i));
        }
        assertEquals(1, DataUtils.peaks.size());
    }

    @Test
    public void processCorrectlyRemovesPeakAfterProcessed(){
        setupPeakTimeSeriesAndAxis();
        int i;
        for(i=0; i<50; i++){
            float[] values = {0, 0, 0};
            DataUtils.process(values, DataUtils.POLLING_RATE*i);
        }
        DataUtils.process(new float[] {100,100,100}, DataUtils.POLLING_RATE*i++);
        int j = i;
        for(; i<j+100; i++){
            float[] values = {0, 0, 0};
            DataUtils.process(values, DataUtils.POLLING_RATE*(i));
        }
        assertEquals(0, DataUtils.peaks.size());
    }

    @Test
    public void testProcessRerunsSGFilterOnNewDataEntersWindow() throws InterruptedException {
        
        setupPeakTimeSeriesAndAxis();
        int i;
        for(i=0; i<50; i++){
            float[] values = {(float) Math.random(), (float) Math.random(), (float) Math.random()};
            DataUtils.process(values, DataUtils.POLLING_RATE*i);
        }
        ArrayList<Float> data = res.processed.get(DataUtils.majorAxisIndex);
        float[] oldValues = {
                data.get(data.size()-1),
                data.get(data.size()-2),
                data.get(data.size()-3),
                data.get(data.size()-4)
        };

        float[] values = {(float) Math.random(), (float) Math.random(), (float) Math.random()};
        DataUtils.process(values, DataUtils.POLLING_RATE*i++);
        DataUtils.process(values, DataUtils.POLLING_RATE*i++);
        DataUtils.process(values, DataUtils.POLLING_RATE*i);

        float[] newValues = {
                data.get(data.size()-4),
                data.get(data.size()-5),
                data.get(data.size()-6),
                data.get(data.size()-7)
        };
        /*
            oldValues[3] and newValues[3] should be equal because its window ends as soon as the
            center node is processed. When that happens is the last time it runs through the filter.
         */
        assertNotEquals(oldValues[0], newValues[0], 0.000001);
        assertNotEquals(oldValues[1], newValues[1], 0.000001);
        assertNotEquals(oldValues[2], newValues[2], 0.000001);
        assertEquals(oldValues[3], newValues[3], 0.000001);
    }

    @Test
    public void testProcessDoesNotRerunSGFilterOnNewDataEntersOutsideWindow() throws InterruptedException {
        
        setupPeakTimeSeriesAndAxis();
        int i;
        for(i=0; i<50; i++){
            float[] values = {(float) Math.random(), (float) Math.random(), (float) Math.random()};
            DataUtils.process(values, DataUtils.POLLING_RATE*i);
        }
        float oldValue = res.processed.get(DataUtils.majorAxisIndex)
                .get(res.processed.get(DataUtils.majorAxisIndex).size()-5);

        float[] values = {(float) Math.random(), (float) Math.random(), (float) Math.random()};
        DataUtils.process(values, DataUtils.POLLING_RATE*i);

        float newValue = res.processed.get(DataUtils.majorAxisIndex)
                .get(res.processed.get(DataUtils.majorAxisIndex).size()-6);
        assertEquals(oldValue, newValue, 0.000001);
    }

    @Test
    public void processDoesNotRemovePeakBeforeProcessed(){
        setupPeakTimeSeriesAndAxis();
        int i;
        for(i=0; i<50; i++){
            float[] values = {0, 0, 0};
            DataUtils.process(values, DataUtils.POLLING_RATE*i);
        }
        DataUtils.process(new float[] {100,100,100}, DataUtils.POLLING_RATE*i++);
        int j = i;
        for(; i<j+30; i++){
            float[] values = {0, 0, 0};
            DataUtils.process(values, DataUtils.POLLING_RATE*i);
        }
        assertEquals(1, DataUtils.peaks.size());
    }

    @Test
    public void testProcessRemovesOverlappingPeaks(){
        setupPeakTimeSeriesAndAxis();
        int i;
        for(i=0; i<50;i++){
            float[] values = {0, 0, 0};
            DataUtils.process(values, DataUtils.POLLING_RATE*i);
        }
        DataUtils.process(new float[] {100,100,100}, DataUtils.POLLING_RATE*i++);
        DataUtils.process(new float[] {0,0,0}, DataUtils.POLLING_RATE*i++);
        DataUtils.process(new float[] {100,100,100}, DataUtils.POLLING_RATE*i++);
        DataUtils.process(new float[] {0,0,0}, DataUtils.POLLING_RATE*i++);
        assertEquals(2, DataUtils.peaks.size());
        int j = i;
        for(; i<j+35;i++){
            float[] values = {0, 0, 0};
            DataUtils.process(values, DataUtils.POLLING_RATE*i);
        }
        assertEquals(0, DataUtils.peaks.size());
    }

    @Test
    public void testProcessDoesNotRemoveNonOverlappingPeaks(){
        setupPeakTimeSeriesAndAxis();
        int i;
        for(i=0; i<50; i++){
            float[] values = {0, 0, 0};
            DataUtils.process(values, DataUtils.POLLING_RATE*i);
        }
        DataUtils.process(new float[] {100,100,100}, DataUtils.POLLING_RATE*i++);

        int j = i;
        for(; i<j+31; i++){
            float[] values = {0, 0, 0};
            DataUtils.process(values, DataUtils.POLLING_RATE*i);
        }
        DataUtils.process(new float[] {100,100,100}, DataUtils.POLLING_RATE*i++);
        assertEquals(2, DataUtils.peaks.size());

        j = i;
        for(; i<j+5; i++){ //TODO: /4->5, need to fix time series data to fix
            float[] values = {0, 0, 0};
            DataUtils.process(values, DataUtils.POLLING_RATE*i);
        }
        assertEquals(1, DataUtils.peaks.size());
    }

    @Test
    public void testProcessIgnoresDataOnNonMajorAxisIndex(){
        setupPeakTimeSeriesAndAxis();
        for(int i=0; i<50; i++){
            float[] values = {(float) Math.random(), (float) Math.random(), (float) Math.random()};
            DataUtils.process(values, DataUtils.POLLING_RATE*i);
        }
        assertEquals(0, res.processed.get(1).size());
    }

    void setupPeakTimeSeriesAndAxis(){
        res.setupPeakTimeSeriesAndAxis();
    }


}
