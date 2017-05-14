package edu.temple.gymminder;

import com.fastdtw.timeseries.TimeSeries;
import com.fastdtw.timeseries.TimeSeriesBase;

import org.junit.*;
import org.junit.Test;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Random;

import edu.temple.gymminder.processing.DataUtils;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

/**
 * Created by rober_000 on 4/12/2017.
 */

public class PeakDetectionTest {

    @Rule
    public final DataUtilsResources res = new DataUtilsResources();

    @Test
    public void testMovingZScorePeakDetectionReturnsOnlyPeakValue() throws NoSuchFieldException, IllegalAccessException {
        Class peakClass = res.getAccessibleDataUtilsClass("Peak");
        for(int i=0; i<50; i++) res.processed.get(0).add(1f);
        res.processed.get(0).add(7f);
        Object o = DataUtils.movingZScorePeakDetection(5, 20, 5, 0);
        float amp = (float) peakClass.getDeclaredField("amplitude").get(o);
        int index = (int) peakClass.getDeclaredField("index").get(o);
        assertEquals(7, amp, 0.0001);
        assertEquals(50, index);
    }

    @Test
    public void testZScorePeakDetectionReturnsOnlyPeakValue() throws NoSuchFieldException, IllegalAccessException {
        Class peakClass = res.getAccessibleDataUtilsClass("Peak");
        TimeSeriesBase.Builder builder = new TimeSeriesBase.Builder();
        for(int i=0; i<50; i++) builder.add(i, 1f);
        builder.add(51, 7f);
        TimeSeries t1 = builder.build();
        ArrayList list = DataUtils.zScorePeakDetection(t1);
        float amp = (float) peakClass.getDeclaredField("amplitude").get(list.get(0));
        int index = (int) peakClass.getDeclaredField("index").get(list.get(0));
        assertEquals(7, amp, 0.0001);
        assertEquals(50, index);
        assertEquals(1, list.size());
    }

    @Test
    public void testMovingZScorePeakDetectionReturnsFinalPeakValue() throws NoSuchFieldException, IllegalAccessException {
        Class peakClass = res.getAccessibleDataUtilsClass("Peak");
        for(int i=0; i<50; i++) res.processed.get(0).add(1f);
        res.processed.get(0).add(7f);
        res.processed.get(0).add(8f);
        Object o = DataUtils.movingZScorePeakDetection(5, 20, 5, 0);
        float amp = (float) peakClass.getDeclaredField("amplitude").get(o);
        int index = (int) peakClass.getDeclaredField("index").get(o);
        assertEquals(8, amp, 0.0001);
        assertEquals(51, index);
    }

    @Test
    public void testZScorePeakDetectionReturnsFinalPeakValue() throws NoSuchFieldException, IllegalAccessException {
        Class peakClass = res.getAccessibleDataUtilsClass("Peak");
        TimeSeriesBase.Builder builder = new TimeSeriesBase.Builder();
        for(int i=0; i<50; i++) builder.add(i, 1f);
        builder.add(51, 7f);
        builder.add(52, 8f);
        TimeSeries t1 = builder.build();
        ArrayList list = DataUtils.zScorePeakDetection(t1);
        float amp = (float) peakClass.getDeclaredField("amplitude").get(list.get(0));
        int index = (int) peakClass.getDeclaredField("index").get(list.get(0));
        assertEquals(8, amp, 0.0001);
        assertEquals(51, index);
        assertEquals(1, list.size());
    }

    @Test
    public void testMovingZScorePeakDetectionDoesNotReturnFinalMaxValueOutsideOfPeak() throws NoSuchFieldException, IllegalAccessException {
        Class peakClass = res.getAccessibleDataUtilsClass("Peak");
        for(int i=0; i<50; i++) res.processed.get(0).add(1f);
        res.processed.get(0).add(7f);
        res.processed.get(0).add(1f);
        res.processed.get(0).add(8f);
        Object o = DataUtils.movingZScorePeakDetection(5, 20, 5, 0);
        float amp = (float) peakClass.getDeclaredField("amplitude").get(o);
        int index = (int) peakClass.getDeclaredField("index").get(o);
        assertEquals(7, amp, 0.0001);
        assertEquals(50, index);
    }

    @Test
    public void testMovingZScorePeakDetectionReturnsNullWithNoPeakValues()  {
        for(int i=0; i<50; i++) res.processed.get(0).add(1f);
        Object o = DataUtils.movingZScorePeakDetection(5, 20, 5, 0);
        assertNull(o);
    }

    @Test
    public void testZScorePeakDetectionReturnsEmptyListWithNoPeakValues(){
        TimeSeriesBase.Builder builder = new TimeSeriesBase.Builder();
        for(int i=0; i<50; i++) builder.add(i, 1f);
        TimeSeries t1 = builder.build();
        ArrayList list = DataUtils.zScorePeakDetection(t1);
        assertEquals(0, list.size());
    }

    @Test
    public void testMovingZScorePeakDetectionDoesNotReturnPeakWithGradualIncrease(){
        for(int i=0; i<50; i++) res.processed.get(0).add(1f);
        Object o = DataUtils.movingZScorePeakDetection(5, 20, 5, 0);
        assertNull(o);
    }

    @Test
    public void testZScorePeakDetectionReturnsListOfSize1WithGradualIncrease(){
        TimeSeriesBase.Builder builder = new TimeSeriesBase.Builder();
        for(int i=0; i<50; i++) builder.add(i, i);
        TimeSeries t1 = builder.build();
        ArrayList list = DataUtils.zScorePeakDetection(t1);
        assertEquals(1, list.size());
    }

    @Test
    public void testMovingZScorePeakDetectionDoesNotUpdateZAfterFoundPeak() throws NoSuchFieldException, IllegalAccessException {
        Class peakClass = res.getAccessibleDataUtilsClass("Peak");
        for(int i=0; i<50; i++) res.processed.get(0).add(1f);
        for(int i=0; i<50; i++) res.processed.get(0).add(50f);
        res.processed.get(0).add(51f);
        Object o = DataUtils.movingZScorePeakDetection(5, 20, 5, 0);

        float amp = (float) peakClass.getDeclaredField("amplitude").get(o);
        int index = (int) peakClass.getDeclaredField("index").get(o);
        assertEquals(51, amp, 0.0001);
        assertEquals(100, index);
    }

    @Test
    public void testZScorePeakDetectionReturnsListWithSizeTwoWithTwoDistinctPeaks(){
        TimeSeriesBase.Builder builder = new TimeSeriesBase.Builder();
        for(int i=0; i<50; i++) builder.add(i, 1f);
        builder.add(51, 7f);
        builder.add(52, 1f);
        builder.add(53, 8f);
        TimeSeries t1 = builder.build();
        ArrayList list = DataUtils.zScorePeakDetection(t1);
        assertEquals(2, list.size());
    }

    @Test
    public void testMovingZScorePeakDetectionReturnsNullWithNegativePeak(){
        int lag=5, window=20, start=50;
        for(int i=0; i<50; i++) res.processed.get(0).add(1f);
        res.processed.get(0).add(-50f);
        for(int i=0; i<50; i++) res.processed.get(0).add(1f);
        Object o = DataUtils.movingZScorePeakDetection(lag, window, 10, 0);
        assertNull(o);
    }

    @Test
    public void testZScorePeakDetectionReturnsEmptyListWithNegativePeak(){
        TimeSeriesBase.Builder builder = new TimeSeriesBase.Builder();
        for(int i=0; i<50; i++) builder.add(i, 1f);
        builder.add(51, -50f);
        TimeSeries t1 = builder.build();
        ArrayList list = DataUtils.zScorePeakDetection(t1);
        assertEquals(0, list.size());
    }

    @Test
    public void testMovingZScorePeakDetectionWorksWithNoiseValues() throws NoSuchFieldException, IllegalAccessException {
        Class peakClass = res.getAccessibleDataUtilsClass("Peak");
        for(int i=0; i<50; i++) res.processed.get(0).add((float) Math.random());
        for(int i=0; i<50; i++) res.processed.get(0).add(50f);
        res.processed.get(0).add(51f);
        Object o = DataUtils.movingZScorePeakDetection(5, 20, 50, 0);
        float amp = (float) peakClass.getDeclaredField("amplitude").get(o);
        int index = (int) peakClass.getDeclaredField("index").get(o);
        assertEquals(51, amp, 0.0001);
        assertEquals(100, index);
    }

    @Test
    public void testZScorePeakDetectionWorksWithNoiseValues() throws NoSuchFieldException, IllegalAccessException {
        Class peakClass = res.getAccessibleDataUtilsClass("Peak");
        TimeSeriesBase.Builder builder = new TimeSeriesBase.Builder();
        for(int i=0; i<50; i++) builder.add(i, Math.random());
        builder.add(51, 7f);
        TimeSeries t1 = builder.build();
        ArrayList list = DataUtils.zScorePeakDetection(t1);
        float amp = (float) peakClass.getDeclaredField("amplitude").get(list.get(0));
        int index = (int) peakClass.getDeclaredField("index").get(list.get(0));
        assertEquals(7, amp, 0.0001);
        assertEquals(50, index);
    }


    @Test
    public void testMovingZScorePeakDetectionReturnsNullWithPeakBeforeStart() {
        int lag=5, window=20, start=50;
        for(int i=0; i<50; i++) res.processed.get(0).add(1f);
        res.processed.get(0).add(50f);
        for(int i=0; i<50; i++) res.processed.get(0).add(1f);
        Object o = DataUtils.movingZScorePeakDetection(lag, window, 10, start+window-lag+1);
        assertNull(o);
    }

    @Test
    public void testMovingZScorePeakDetectionReturnsObjectWithPeakAtStartOfStart(){
        int lag=5, window=20, start=50;
        for(int i=0; i<50; i++) res.processed.get(0).add(1f);
        res.processed.get(0).add(50f);
        for(int i=0; i<50; i++) res.processed.get(0).add(1f);
        Object o = DataUtils.movingZScorePeakDetection(lag, window, 10, start+window-lag);
        assertNotNull(o);
    }

    @Test
    public void testZScorePeakDetectionReturnsPeakAtStartOfTimeSeries() throws NoSuchFieldException, IllegalAccessException {
        Class peakClass = res.getAccessibleDataUtilsClass("Peak");
        TimeSeriesBase.Builder builder = new TimeSeriesBase.Builder();
        builder.add(0, 7f);
        for(int i=1; i<=50; i++) builder.add(i, 1f);
        TimeSeries t1 = builder.build();
        ArrayList list = DataUtils.zScorePeakDetection(t1);
        float amp = (float) peakClass.getDeclaredField("amplitude").get(list.get(0));
        int index = (int) peakClass.getDeclaredField("index").get(list.get(0));
        assertEquals(7, amp, 0.0001);
        assertEquals(0, index);
        assertEquals(1, list.size());
    }

    @Test
    public void testZScorePeakDetectionReturnsPeakInMiddleOfTimeSeries() throws NoSuchFieldException, IllegalAccessException {
        Class peakClass = res.getAccessibleDataUtilsClass("Peak");
        TimeSeriesBase.Builder builder = new TimeSeriesBase.Builder();
        for(int i=0; i<=22; i++) builder.add(i, 1f);
        builder.add(23, 7f);
        for(int i=24;i<50;i++) builder.add(i, 1f);
        TimeSeries t1 = builder.build();
        ArrayList list = DataUtils.zScorePeakDetection(t1);
        float amp = (float) peakClass.getDeclaredField("amplitude").get(list.get(0));
        int index = (int) peakClass.getDeclaredField("index").get(list.get(0));
        assertEquals(7, amp, 0.0001);
        assertEquals(23, index);
        assertEquals(1, list.size());
    }

    @Test
    public void testSuccsesfulRemovalOfPeaks() throws NoSuchMethodException, IllegalAccessException, InvocationTargetException, InstantiationException, NoSuchFieldException {
        Class peakClass = res.getAccessibleDataUtilsClass("Peak");
        Constructor constr = peakClass.getDeclaredConstructors()[0];
        constr.setAccessible(true);
        Object originalPeak = constr.newInstance(0, 10f);
        ArrayList<Object> peaks = new ArrayList<>();
        Random rand = new Random();
        for(int i=0;i<10;i++){
            peaks.add(constr.newInstance(0, rand.nextFloat()*3));
        }
        assertEquals(10, peaks.size());
        peaks = (ArrayList<Object>) DataUtils.class
                .getDeclaredMethod("reducePeaks", ArrayList.class, originalPeak.getClass())
                .invoke(null, peaks, originalPeak);
        assertEquals(0, peaks.size());
    }

    @Test
    public void testNonRemovalOfValidPeaks() throws IllegalAccessException, InvocationTargetException, InstantiationException, NoSuchMethodException, NoSuchFieldException {
        Class peakClass = res.getAccessibleDataUtilsClass("Peak");
        Constructor constr = peakClass.getDeclaredConstructors()[0];
        constr.setAccessible(true);
        Object originalPeak = constr.newInstance(0, 10f);
        ArrayList<Object> peaks = new ArrayList<>();
        Random rand = new Random();
        for(int i=0;i<10;i++){
            peaks.add(constr.newInstance(0, 6.68f+rand.nextFloat()*3f));
        }
        assertEquals(10, peaks.size());
        peaks = (ArrayList<Object>) DataUtils.class
                .getDeclaredMethod("reducePeaks", ArrayList.class, originalPeak.getClass())
                .invoke(null, peaks, originalPeak);
        assertEquals(10, peaks.size());
    }


}
