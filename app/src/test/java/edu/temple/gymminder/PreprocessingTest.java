package edu.temple.gymminder;

import org.junit.Rule;
import org.junit.Test;

import java.util.ArrayList;

import edu.temple.gymminder.processing.DataUtils;

import static org.junit.Assert.assertEquals;

/**
 * Created by rober_000 on 4/12/2017.
 */

public class PreprocessingTest {

    @Rule
    public final DataUtilsResources res = new DataUtilsResources();

    @Test
    public void applySavitzkyGolayFilterDoesNotCrash() {
        ArrayList<Float> data = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            data.add(1f);
        }
        DataUtils.applySavitzkyGolayFilter(data);
    }

    @Test
    public void realTimeSavitzkyGolayFilterDoesNotCrash() {
        for (int i = 0; i < 10; i++) {
            res.data.get(0).add((float) Math.random());
        }
        for (int i = 0; i < 10; i++) {
            DataUtils.applySGFilterRealtime(i, res.data.get(0), res.processed.get(0));
        }
    }

    @Test
    public void testSavitzkyGolayFilter() {
        ArrayList<Float> data = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            data.add((float) i);
        }
        data = DataUtils.applySavitzkyGolayFilter(data);
        float[] expected = {.28f, .95f, 1.9f, 3f, 4f, 5f, 6f, 7.95f, 7.61f, 5.71f};
        for (int i = 0; i < 10; i++) {
            assertEquals(expected[i], data.get(i), .01);
        }
    }

    @Test
    public void testSavitzkyGolayFilterRealtime() {
        ArrayList<Float> data = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            data.add((float) i);
        }
        for (int i = 0; i < 10; i++) {
            DataUtils.applySGFilterRealtime(i, data, res.processed.get(0));
        }
        float[] expected = {.28f, .95f, 1.9f, 3f, 4f, 5f, 6f, 7.95f, 7.61f, 5.71f};
        for (int i = 0; i < 10; i++) {
            assertEquals(expected[i], res.processed.get(0).get(i), .01);
        }
        //Run filter again to simulate re-processing data points at each run of DataUtils.process()
        for (int i = 0; i < 10; i++) {
            DataUtils.applySGFilterRealtime(i, data, res.processed.get(0));
        }
        for (int i = 0; i < 10; i++) {
            assertEquals(expected[i], res.processed.get(0).get(i), .01);
        }
    }

    @Test
    public void testInterpolation(){
        ArrayList<ArrayList<Float>> data = new ArrayList<>();
        data.add(new ArrayList<Float>());
        data.get(0).add(10f);
        //reinitialize DataUtils to use our new data sets
        DataUtils.init(data, null, null);
        float x = DataUtils.interpolate(20, .2f, 0);
        assertEquals(15, x, .01);
        x = DataUtils.interpolate(30, .2f, 0);
        assertEquals(20, x, .01);
        x = DataUtils.interpolate(30, .3f, 0);
        assertEquals(16.66, x, .01);
    }

    @Test
    public void testAverage(){
        float[] avgNode = {10f, 10f};
        float[] newNode = {20f, 10f};
        float[] result = DataUtils.average(avgNode, newNode[0], newNode[1]);
        assertEquals(15, result[0], .1);
        assertEquals(20, result[1], .1);
        newNode = new float[]{20f, 20f};
        avgNode = new float[]{10f, 10f};
        result = DataUtils.average(avgNode, newNode[0], newNode[1]);
        assertEquals(16.66f, result[0], .1);
        assertEquals(30f, result[1], .1);
    }

}
