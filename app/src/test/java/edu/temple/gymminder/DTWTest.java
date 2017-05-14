package edu.temple.gymminder;

import com.fastdtw.dtw.FastDTW;
import com.fastdtw.dtw.WarpPath;
import com.fastdtw.timeseries.TimeSeries;
import com.fastdtw.timeseries.TimeSeriesBase;
import com.fastdtw.util.Distances;

import org.junit.Rule;
import org.junit.Test;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

import edu.temple.gymminder.processing.DataUtils;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

/**
 * Created by rober_000 on 4/12/2017.
 */

public class DTWTest {

    @Rule
    public final DataUtilsResources res = new DataUtilsResources();

    @Test
    public void testGetLastMatchingIndexOfFirst(){
        TimeSeriesBase.Builder builder1 = TimeSeriesBase.builder();
        TimeSeriesBase.Builder builder2 = TimeSeriesBase.builder();
        for(int i=0; i<2; i++){
            builder1 = builder1.add(i, 0);
        }
        for(int i=0; i<10; i++){
            builder1 = builder1.add(i+2, i);
            builder2 = builder2.add(i, i);
        }
        TimeSeries t1 = builder1.build();
        TimeSeries t2 = builder2.build();
        WarpPath path = FastDTW.compare(t1, t2, Distances.EUCLIDEAN_DISTANCE).getPath();
        int index = DataUtils.getLastMatchingIndexOfFirst(path);
        //Last matching index should be third, since we added 2 extra values
        assertEquals(2, index);
    }

    @Test
    public void testGetFirstMatchingIndexOfLast(){
        TimeSeriesBase.Builder builder1 = TimeSeriesBase.builder();
        TimeSeriesBase.Builder builder2 = TimeSeriesBase.builder();
        for(int i=0; i<10; i++){
            builder1 = builder1.add(i, i);
            builder2 = builder2.add(i, i);
        }
        for(int i=0; i<2; i++){
            builder1 = builder1.add(i+10, 10);
        }
        TimeSeries t1 = builder1.build();
        TimeSeries t2 = builder2.build();
        WarpPath path = FastDTW.compare(t1, t2, Distances.EUCLIDEAN_DISTANCE).getPath();
        int index = DataUtils.getFirstMatchingIndexOfLast(path);
        assertEquals(9, index);
    }

    @Test
    public void testSubSeries(){
        TimeSeriesBase.Builder builder = TimeSeriesBase.builder();
        for(int i=0; i<10; i++){
            builder = builder.add(i, i);
        }
        TimeSeries ts = builder.build();
        assertEquals(10, ts.size());
        ts = DataUtils.subSeries(ts, 2, 8);
        assertEquals(2, ts.getMeasurement(0, 0), 0);
        assertEquals(7, ts.getMeasurement(5, 0), 0);
        assertEquals(6, ts.size());
    }

    @Test
    public void testCalcFeatures() throws IllegalAccessException, InvocationTargetException, InstantiationException, NoSuchMethodException, NoSuchFieldException {
        Class peakClass = res.getAccessibleDataUtilsClass("Peak");
        TimeSeriesBase.Builder builder = TimeSeriesBase.builder();
        for (int i = 0; i < 10; i++) {
            builder = builder.add(i, i);
        }
        TimeSeries t1 = builder.build();
        TimeSeries t2 = builder.build();
        Constructor constr = peakClass.getDeclaredConstructors()[0];
        constr.setAccessible(true);
        Object originalPeak = constr.newInstance(0, 9f);
        double[] features = (double[]) DataUtils.class
                .getMethod("calcFeatures", TimeSeries.class, originalPeak.getClass(), TimeSeries.class)
                .invoke(null, t1, originalPeak, t2);
        double[] expected = {0, 9, 0, 2.87228, 5.33853};
        assertArrayEquals(expected, features, .001);
    }

    @Test
    public void testAccept() throws IllegalAccessException, InvocationTargetException, InstantiationException, NoSuchMethodException, NoSuchFieldException {
        Class boundsClass = res.getAccessibleDataUtilsClass("DetectedBounds");
        Constructor constr = boundsClass.getDeclaredConstructors()[0];
        constr.setAccessible(true);
        //TODO: change constructor call when coefficients are found
        Object bounds = constr.newInstance(0, 0, 0f, 0f, 0f, 0f, 0f);
        boolean x = (boolean) DataUtils.class
                .getDeclaredMethod("accept", bounds.getClass())
                .invoke(null, bounds);
        assertEquals(true, x);

        bounds = constr.newInstance(-1, 0, 0f, 0f, 0f, 0f, 0f);
        x = (boolean) DataUtils.class
                .getDeclaredMethod("accept", bounds.getClass())
                .invoke(null, bounds);
        assertEquals(x, true);
    }

}
