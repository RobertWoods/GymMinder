package edu.temple.gymminder;

import org.junit.Rule;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.mockito.Mock;
import org.mockito.runners.MockitoJUnitRunner;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;

import edu.temple.gymminder.processing.DataUtils;

import static org.junit.Assert.*;
import static org.mockito.Mockito.when;

/**
 * Example local unit test, which will execute on the development machine (host).
 *
 * @see <a href="http://d.android.com/tools/testing">Testing documentation</a>
 */
@RunWith(MockitoJUnitRunner.class)
public class BasicDataUtilsMethodsTest {

    @Rule
    public final DataUtilsResources res = new DataUtilsResources();

    @Test
    public void testDataUtilsConstantRiemann() {
        ArrayList<Float> data = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            data.add(1f);
        }
        float f[] = DataUtils.riemann(data);
        for (int i = 0; i < 9; i++) {
            assertEquals(1, f[i], .1);
        }
    }

    @Test
    public void testDataUtilsConstantIntegration() {
        ArrayList<Float> data = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            data.add(1f);
        }
        float f = DataUtils.sum(DataUtils.riemann(data));
        assertEquals(9, f, 1);
    }

    @Test
    public void testDataUtilsLinearRiemann() {
        ArrayList<Float> data = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            data.add(1f * i);
        }
        float f[] = DataUtils.riemann(data);
        for (int i = 0; i < 9; i++) {
            assertEquals(i, f[i], .1);
        }
    }

    @Test
    public void testDataUtilsLinearIntegration() {
        ArrayList<Float> data = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            data.add(1f * i);
        }
        float f = DataUtils.sum(DataUtils.riemann(data));
        assertEquals(36, f, 1);
    }

    @Test
    public void testSumWithEndIndex(){
        float floats[] = new float[20];
        for(int i=0; i<20; i++){
            floats[i] = i;
        }
        int sum = 0;
        for(int i=0; i<20; i++){
            assertEquals(sum, DataUtils.sum(floats, i), 0.00001);
            sum+=i;
        }
    }

    @Test
    public void testPartialSums(){
        float floats[] = new float[20];
        for(int i=0; i<20; i++){
            floats[i] = i;
        }
        float sums[] = {0, 1, 3, 6,10,15,21,28,36,45, 55,66,78,91,105,120,136,153,171,190};
        float partials[] = DataUtils.partialSums(floats);
        assertArrayEquals(sums, partials, 0.000001f);
    }

    @Mock
    BufferedReader reader;

    @Test
    public void testLoadRepetitionPatternTimeSeries() throws IOException, NoSuchFieldException, IllegalAccessException {
        when(reader.readLine())
                .thenReturn("1.4,2,3.6,4,5,6.3,7,8,9")
                .thenReturn("5,6.3")
                .thenReturn("0");
        DataUtils.loadRepetitionPatternTimeSeries(reader);
        //TODO: use reflection and set these back to private... do it with some other things too...
        assertNotNull(DataUtils.repPeak);
        assertNotNull(DataUtils.repTimeSeries);
        assertEquals(9, DataUtils.repTimeSeries.size());
        Object o = DataUtils.repPeak;
        assertEquals(5, res.getAccessibleDataUtilsClass("Peak").getDeclaredField("index").get(o));
        assertEquals(6.3f, (float) res.getAccessibleDataUtilsClass("Peak").getDeclaredField("amplitude").get(o), 0.01);
        assertEquals(0, DataUtils.majorAxisIndex);
    }

}
    

