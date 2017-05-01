package edu.temple.gymminder;

import android.content.Context;
import android.util.Log;
import android.util.SparseArray;

import com.fastdtw.dtw.FastDTW;
import com.fastdtw.dtw.TimeWarpInfo;
import com.fastdtw.dtw.WarpPath;
import com.fastdtw.matrix.ColMajorCell;
import com.fastdtw.timeseries.TimeSeries;
import com.fastdtw.timeseries.TimeSeriesBase;
import com.fastdtw.util.Distances;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import edu.temple.gymminder.models.Repetition;

/**
 * Created by rober_000 on 2/10/2017.
 */

public class DataUtils {

    //TODO: Refactor at least data processing part to instantiated class

    public static final float MS2S_CONVERSION = 1.0f / 1000000.0f;
    public static final long SECOND = 1000000;
    public static final float[] SG_FILTER = {-2, 3, 6, 7, 6, 3, -2};
    public static final float FILTER_SUM = sum(SG_FILTER);
    private static final float PERIOD = .1f;
    private static final double EXPANSION_VALUE = 1.5;
    public static final long POLLING_FREQUENCY = 30;
    public static final long POLLING_RATE = SECOND / POLLING_FREQUENCY;
    private static final double PEAK_SIMILARITY_FACTOR = 1.5;
    private static final long ERROR = 1000;

    private static float[] avgNode = null;
    private static ArrayList<ArrayList<Float>> data;
    private static ArrayList<ArrayList<Float>> processedData;
    private static ArrayList<Long> timestamps;
    public static SparseArray<Peak> peaks = new SparseArray<>();

    public static Peak repPeak;
    public static TimeSeries repTimeSeries;
    public static int majorAxisIndex;

    private static ExecutorService executorService;

    private static Listener listener;

    static void init(ArrayList<ArrayList<Float>> dataList, ArrayList<Long> time) {
        data = dataList;
        timestamps = time;
        processedData = new ArrayList<>(3);
        peaks = new SparseArray<>();
        for (int i = 0; i < 3; i++) processedData.add(new ArrayList<Float>());
    }

    static void init(ArrayList<ArrayList<Float>> dataList, ArrayList<Long> time, ArrayList<ArrayList<Float>> processedDataList) {
        init(dataList, time);
        processedData = processedDataList;
        avgNode = null;
    }

    static void init(ArrayList<ArrayList<Float>> dataList, ArrayList<Long> time, ArrayList<ArrayList<Float>> processedDataList, File file) {
        init(dataList, time, processedDataList);
        try {
            loadRepetitionPatternTimeSeries(new BufferedReader(new FileReader(file)));
            peaks = new SparseArray<>((int) (EXPANSION_VALUE * (repTimeSeries.size() - repPeak.index)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    static void setListener(Listener l) {
        listener = l;
        executorService = Executors.newSingleThreadExecutor();
    }

    //Need this to prevent possible memory leak
    static void removeListener() {
        executorService.shutdownNow();
        listener = null;
    }

    /**
     * @param x     x acceleration to be added into array
     * @param y     y acceleration to be added into array
     * @param z     z acceleration to be added into array
     * @param alpha low pass filter parameter
     */
    static void addWithLowPassFilter(float x, float y, float z, float alpha) {
        if (data.get(0).size() > 0) {
            float ox = data.get(0).get(data.get(0).size() - 1);
            float oy = data.get(1).get(data.get(0).size() - 1);
            float oz = data.get(2).get(data.get(0).size() - 1);
            data.get(0).add(ox + alpha * (x - ox));
            data.get(1).add(oy + alpha * (y - oy));
            data.get(2).add(oz + alpha * (z - oz));
        } else {
            addUnfiltered(x, y, z);
        }
    }

    static void addUnfiltered(float x, float y, float z) {
        data.get(0).add(x);
        data.get(1).add(y);
        data.get(2).add(z);
    }

    /**
     * This is an intermediate step in approximating the definite integral of data points, for
     * use with sum
     *
     * @param list list of float values for which to calculate riemann rectangles
     * @return array of riemann rectangles
     */
    static float[] riemann(List<Float> list) {
        float[] velocity = new float[list.size()-1];
        Iterator<Float> iterator = list.listIterator();
        int i = 0;
        while (iterator.hasNext()) {
            float value = iterator.next();
            if(!iterator.hasNext()) return velocity;
            velocity[i] = value;
            i++;
        }
        return velocity;
    }

    /**
     * @param floats list of floats to be summed
     * @param end    index at which to stop summation
     * @return the sum of values in floats
     */
    static float sum(float[] floats, int end) {
        float sum = 0f;
        for(int i=0; i<end;i++){
            sum+=floats[i];
        }
        return sum;
    }

    static float sum(float[] floats) {
        return sum(floats, floats.length);
    }

    /**
     *
     * @param floats    An array of values representing riemann sums
     * @return          List of partial sums at each value
     */
    static float[] partialSums(float[] floats){
        float[] data = new float[floats.length];
        for(int i=0; i<floats.length; i++){
            data[i] = (sum(floats, i+1));
        }
        return data;
    }

    /**
     * @param floats array of floats for which to find the maximum and average value
     * @return a float array whose first value is the max and whose second value is the average
     */
    static float[] maxAndAvg(float[] floats) {
        float f[] = new float[2];
        float max = floats[0];
        float sum = 0;
        for (int i = 0; i < floats.length; i++) {
            max = max > floats[i] ? max : floats[i];
            sum += floats[i];
        }
        f[0] = max;
        f[1] = sum / floats.length;
        return f;
    }

    /**
     * Adds data point into the uniformly-spaced and smoothed data array, and processes peak values
     * <p>
     * avgNode is used by this method when the real polling frequency of the device
     * is higher than the desired polling frequency, and represents all data points not yet
     * placed into the data array because the time difference is still too small where
     * avgNode[0] = value, and avgNode[1] = duration
     * <p>
     * If the polling frequency of the device is lower than the desired polling frequency it
     * interpolates a data point at the desired frequency using the newest data point and the
     * data point most recently added into the data array
     *
     * @param values    event values
     * @param timestamp timestamp of event
     */
    static void process(final float[] values, final long timestamp) {
        executorService.submit(new Runnable() {
            @Override
            public void run() {
                int i = majorAxisIndex;

                float x = Math.abs(values[i]) > 0.015 ? values[i] : 0;
                //First time adding a node, just add it lel
                if (timestamps.size() == 0) {
                    timestamps.add(timestamp);
                    data.get(i).add(x);
                    processedData.get(i).add(x);
                    return;
                }
                float duration = (timestamp - timestamps.get(timestamps.size() - 1)) * MS2S_CONVERSION;
                long longDuration = timestamp - timestamps.get(timestamps.size()-1);

                if ((longDuration + ERROR) < POLLING_RATE || avgNode != null) {
                    //average the points with sum node
                    avgNode = average(avgNode, x, duration);
                    duration = avgNode[1];
                    x = avgNode[0];
                }
                if ((longDuration + ERROR) >= PERIOD) {
                    //We can approximate timestamp value by adding .1s to previous value
                    //Maybe not the best idea since it (maybe) causes drift when we interpolate, idk :d
                    timestamps.add(timestamps.get(timestamps.size() - 1) + POLLING_RATE);
                    data.get(i).add(x);
                    int size = processedData.get(i).size();
                    for (int j = size; j >= size - (SG_FILTER.length / 2 ) && j > 0; j--) {
                /*
                    We want to re-process any data points that didn't have enough data to the right
                    for the entire filter to run on. The alternative to this is to delay the signal
                    by waiting until we have enough data points for the entire filter.
                 */
                        applySGFilterRealtime(j, data.get(i), processedData.get(i));
                    }
                    avgNode = null;
                    Peak newPeak = detectPeak(size, 1);
                    //TODO: Probably want to put this in a thread that enqueues new peaks to check
                    if (newPeak != null) {
                    /*
                        The key for our HashMap is the index of the new peak, which is the current
                        size of processedData at the time of insertion, plus the difference between
                        the size of repTimeSeries and the index of repPeak multiplied by
                        EXPANSION_VALUE.
                     */
                        peaks.put((int) (newPeak.index + EXPANSION_VALUE *
                                (repTimeSeries.size() - repPeak.index)), newPeak);
                    }
                    if (peaks.get(processedData.get(i).size()) != null) {
                    /*
                        Because we used an index for our HashMap key value, we can use the
                        current index to detect if any peaks are ready to be examined
                     */
                        TimeSeriesBase.Builder builder = TimeSeriesBase.builder();
                        for (int j = 0; j < processedData.get(i).size(); j++)
                            builder = builder.add(j, processedData.get(i).get(j));
                        TimeSeries t1 = builder.build();
                        DetectedBounds bounds = detectBounds(t1, peaks.get(processedData.get(i).size()));
                        int peakIndex = peaks.get(processedData.get(i).size()).index;
                        peaks.remove(processedData.get(i).size());
                        if (accept(bounds)) {
                        /*
                            This was a valid repetition, so we want to vibrate and remove any
                            potential peaks that we now know are contained within the repetition
                         */
                            if (listener != null) listener.respondToRep();
                            for (int j = processedData.get(i).size() + 1; j < (processedData.get(i).size() + 1 + (bounds.e - peakIndex)); j++) {
                                if (peaks.get(j) != null) peaks.remove(j);
                            }
                        }

                    }

                }
            }
        });

    }

    private static Peak detectPeak(int index, Object... args) {
        if (args.length == 0) {
            return detectPeakQRSMethod();
        } else {
            return movingZScorePeakDetection(10, 20, 1.5, processedData.get(majorAxisIndex).size());
        }
    }

    private static Peak detectPeakQRSMethod() {
        //TODO: Implement selection of peak candidates: https://www.ncbi.nlm.nih.gov/pubmed/18269982
        return null;
    }

    /**
     * @param axes An array of event values for each axis
     * @return The index of the major axis in the axes array
     */
    public static int detectMajorAxis(ArrayList<ArrayList<Float>> axes) {
        float maxDifference = 0;
        int index = 0;
        
        // Check the difference between the min and max values for each axis
        // TODO: Add reference reps and change algorithm to use min DTW distance to reference
        int i = 0;
        for (ArrayList<Float> axis : axes) {
            float min = 0, max = 0;
            for (float value : axis) {
                if (value < min)
                    min = value;
                if (value > max)
                    max = value;
            }

            float difference = Math.abs(max - min);
            if (difference > maxDifference) {
                maxDifference = difference;
                index = i;
            }
            i++;
        }
        return index;
    }

    /**
     * @param lag    The amount of data points to use to calculate std and mean
     * @param window The total number of data points we want to check
     * @param z      The z-score threshold for detection
     * @param start  The beginning of the detection window
     * @return A new peak representing the max peak in the data set, or null if none
     */
    public static Peak movingZScorePeakDetection(int lag, int window, double z, int start) {
        //Implementation of: http://stackoverflow.com/q/22583391/
        //TODO: Because we filter our data it might be better to run this on the unprocessed data
        //'influence' is 0
        //For now we only detect the peak within the lag
        int m = majorAxisIndex;
        if (processedData.get(m).size() < lag) return null;
        start -= window;
        start = start >= 0 ? start : 0;
        //Calculate std and mean for first lag samples
        double mean = 0;
        for (int i = start; i < start + lag; i++) {
            mean += processedData.get(m).get(i);
        }
        mean /= lag;
        double std = 0;
        for (int i = start; i < start + lag; i++) {
            std += Math.pow(processedData.get(m).get(i) - mean, 2);
        }
        std /= lag;

        //Begin search
        double max = 0;
        int index = -1;
        for (int i = start + lag; i < processedData.get(m).size(); i++) {
            float repPeakValue = repPeak != null ? repPeak.amplitude : 0;
            if (((processedData.get(m).get(i) - mean) / std) > z &&
                    processedData.get(m).get(i) * PEAK_SIMILARITY_FACTOR > repPeakValue) {
                index = processedData.get(m).get(i) > max ? i : index;
                max = processedData.get(m).get(i) > max ? processedData.get(m).get(i) : max;
            } else if (index > 0) {
                //If we have an index but the new data point is not in a peak, we exit
                break;
            } else {
                //Only recalculate mean and std if we are not in a peak
                mean -= (processedData.get(m).get(i - lag) / lag);
                mean += (processedData.get(m).get(i) / lag);
                std = 0;
                for (int j = i - lag + 1; j <= i; j++) {
                    std += Math.pow(processedData.get(m).get(j), 2);
                }
                std /= lag;
            }
        }
        if (index == -1) {
            return null;
        }
        return new Peak(index, (float) max);
    }

    /**
     * @param x        data value to be interpolated
     * @param duration duration > PERIOD since last added node
     * @param i        index of data array being used for interpolation
     * @return interpolated data value
     */
    static float interpolate(float x, float duration, int i) {
        float old = data.get(i).get(data.get(i).size() - 1);
        return old + (PERIOD) * (x - old) / (duration);
    }

    /**
     * @param avgNode     node to be modified
     * @param newValue    new value to be added to average
     * @param newDuration new duration to be added to average
     * @return the modified avgNode
     */
    static float[] average(float[] avgNode, float newValue, float newDuration) {
        if (avgNode == null) {
            avgNode = new float[]{newValue, newDuration};
        } else {
            avgNode[0] = avgNode[0] * (avgNode[1] / (avgNode[1] + newDuration))
                    + newValue * (newDuration / (avgNode[1] + newDuration));
            avgNode[1] = avgNode[1] + newDuration;
        }
        return avgNode;
    }


    /**
     * @param data array of data points to be smoothed
     * @return smoothed data points
     */
    static ArrayList<Float> applySavitzkyGolayFilter(ArrayList<Float> data) {
        ArrayList<Float> filtered = new ArrayList<>(data.size());
        for (int i = 0; i < data.size(); i++) {
            float sum = 0;
            for (int j = 0; j < SG_FILTER.length; j++) {
                if (i + j - SG_FILTER.length / 2 >= 0 && i + j - SG_FILTER.length / 2 < data.size()) {
                    sum += SG_FILTER[j] * data.get(i + j - SG_FILTER.length / 2);
                }
            }
            filtered.add(i, sum / FILTER_SUM);
        }
        return filtered;
    }

    /**
     * @param index         the index of the data point to apply the filter to
     * @param data          data array to be filtered
     * @param processedData data array to be inserted into
     */
    static void applySGFilterRealtime(int index, ArrayList<Float> data, ArrayList<Float> processedData) {
        float sum = 0;
        for (int i = 0; i < SG_FILTER.length; i++) {
            if (i + index - SG_FILTER.length / 2 >= 0 && i + index - SG_FILTER.length / 2 < data.size()) {
                sum += SG_FILTER[i] * data.get(i + index - SG_FILTER.length / 2);
            }
        }
        if (index == processedData.size()) processedData.add(index, sum / FILTER_SUM);
        else processedData.set(index, sum / FILTER_SUM);
    }

    /**
     *
     * @param peaks         List of peaks for which we want to reduce any potential bad peaks
     * @param originalPeak  Original peak, which we used to calculate minimum amplitude threshold
     *                      of other peaks
     * @return A List of peaks that is a subset of the original, with peaks with too small a
     * amplitude removed.
     */
    public static ArrayList<Peak> reducePeaks(ArrayList<Peak> peaks, Peak originalPeak) {
        //TODO: actually call this when detecting peaks
        ArrayList<Peak> result = new ArrayList<>();
        for (Peak p : peaks) {
            if (p.amplitude >= (Math.pow(PEAK_SIMILARITY_FACTOR, -1) * originalPeak.amplitude)) {
                result.add(p);
            }
        }
        return result;
    }

    /**
     *
     * @param series    TimeSeries for which we want to return a subseries
     * @param start     The index at which we want our subseries to begin (inclusive)
     * @param end       The index at which we want our subseries to end (exclusive)
     * @return          A TimeSeries object representing a subsection of series
     */
    public static TimeSeries subSeries(TimeSeries series, int start, int end) {
        TimeSeriesBase.Builder builder = TimeSeriesBase.builder();
        start = start < 0 ? 0 : start;
        for (int i = start; i < end && i < series.size(); i++) {
            builder.add(series.getTimeAtNthPoint(i), series.getMeasurement(i, 0));
        }
        return builder.build();
    }

    public static TimeSeries seriesFromList(List<Float> data){
        TimeSeriesBase.Builder builder = new TimeSeriesBase.Builder();
        int i = 0;
        for(float f : data){
            builder = builder.add(i++, f);
        }
        return builder.build();
    }

    /**
     * This and getFirstMatchingIndexOfLast are used to calculate the bounds of the repetition. This
     * method in is used to calculate the beginning of the rep.
     * @param path path representing the DTW path between two TimeSeries.
     * @return The last of our acceleration stream corresponding to the first index of our
     * reference TimeSeries.
     */
    public static int getLastMatchingIndexOfFirst(WarpPath path) {
        ColMajorCell cell = path.get(0);
        int startIndexI = cell.getRow(); //This might need to be cell.getCol()
        int i;
        for (i = 1; startIndexI == cell.getRow() && i < path.size(); i++) {
            cell = path.get(i);
        }
        //Subtract 2 because the index is not valid, and then i++ was called again
        return path.get(i - 2).getCol(); //If cell.getRow() changes, so does this
    }


    /**
     * This and getFirstMatchingIndexOfLast are used to calculate the bounds of the repetition. This
     * method in is used to calculate the end of the rep.
     * @param path path representing the DTW path between two TimeSeries.
     * @return The first of our acceleration stream corresponding to the last index of our
     * reference TimeSeries.
     */
    public static int getFirstMatchingIndexOfLast(WarpPath path) {
        ColMajorCell cell = path.get(path.size() - 1);
        int endIndexI = cell.getRow(); //ditto comments above
        int i;
        for (i = path.size() - 2; endIndexI == cell.getRow() && i >= 0; i--) {
            cell = path.get(i);
        }
        //Work backwards this time so add the two back
        return path.get(i + 2).getCol();
    }

    /**
     * @param t1     acceleration stream
     * @param t1Peak peak candidate for acceleration stream
     * @return DetectedBounds representing bounds of candidate for repetition
     */
    public static DetectedBounds detectBounds(TimeSeries t1, Peak t1Peak) {
        int s = (int) (t1Peak.index - EXPANSION_VALUE * repPeak.index);
        int e = (int) (t1Peak.index + EXPANSION_VALUE * (repTimeSeries.size() - repPeak.index));
        int s1 = s;
        t1 = subSeries(t1, s, e);
        TimeWarpInfo info = FastDTW.compare(t1, repTimeSeries, Distances.EUCLIDEAN_DISTANCE);
        //Last element s in R -> C[0]
        s = getLastMatchingIndexOfFirst(info.getPath());
        //First element e in R -> C[n]
        e = getFirstMatchingIndexOfLast(info.getPath());
        //Find min, mean, std, rms, dur in R[s':e']
        t1 = subSeries(t1, s, e);
        double[] f = calcFeatures(t1, t1Peak);
        //Add s1 back to start and end index to get absolute start and end indices instead of relative
        return new DetectedBounds(s+s1, e+s1, f[0], f[1], f[2], f[3], f[4]);
    }

    /**
     * @param t1     TimeSeries for which to calculate features
     * @param t1Peak Peak of repetition
     * @param t2     Repetition pattern TimeSeries for which to calculate distance from
     * @return feature values { distance, max, min, standard deviation, root mean square }
     */
    public static double[] calcFeatures(TimeSeries t1, Peak t1Peak, TimeSeries t2) {
        //TODO maybe find a way to get dst faster
        double dst = FastDTW.compare(t1, t2, Distances.EUCLIDEAN_DISTANCE).getDistance();
        double max = t1Peak.amplitude;
        double min = Double.MAX_VALUE;
        double mean = 0;
        double std = 0;
        double rms = 0;
        for (int i = 0; i < t1.size(); i++) {
            double value = t1.getMeasurement(i, 0);
            min = min < value ? min : value;
            mean += value;
            rms += value * value;
        }
        mean /= t1.size();
        for (int i = 0; i < t1.size(); i++) {
            double value = t1.getMeasurement(i, 0);
            std += Math.pow((value - mean), 2);
        }
        rms = Math.sqrt((rms / t1.size()));
        //Using population std I guess o3o
        std = Math.sqrt((std / (t1.size())));
        return new double[]{dst, max, min, std, rms};
    }

    public static ArrayList<Repetition> calculateReps(TimeSeries t1) {
        ArrayList<Peak> peaks = zScorePeakDetection(t1);
        if (repPeak != null) peaks = reducePeaks(peaks, repPeak);
        ArrayList<Repetition> repetitions = new ArrayList<>(peaks.size());
        ArrayList<DetectedBounds> finalBounds = new ArrayList<>(repetitions.size());
        for (Peak p : peaks) {
            DetectedBounds bounds = detectBounds(t1, p);
            if (accept(bounds) && notOverlapping(p, finalBounds)) {
                repetitions.add(repetitionFromPeakAndBounds(p, bounds, t1));
                finalBounds.add(bounds);

            }
        }
        return repetitions;
    }

    private static boolean notOverlapping(Peak peak, ArrayList<DetectedBounds> bounds) {
        for(DetectedBounds b : bounds){
            if(peak.index <= b.e && peak.index >= b.s) return false;
        }
        return true;
    }


    private static Repetition repetitionFromPeakAndBounds(Peak peak, DetectedBounds bounds, TimeSeries t1){
        if(peak==null || bounds==null) return null;
        return new Repetition(bounds.s, bounds.e, peak.index, peak.amplitude, t1);
    }

    public static ArrayList<Peak> zScorePeakDetection(TimeSeries t1) {
        double z = 1;
        double sum = 0;
        for(int i=0;i<t1.size();i++){
            sum+=t1.getMeasurement(i, 0);
        }
        double mean = sum / t1.size();
        sum = 0;
        for(int i=0;i<t1.size();i++){
            sum += Math.pow(t1.getMeasurement(i, 0) - mean, 2);
        }
        double std = Math.sqrt(sum / t1.size());
        ArrayList<Integer> indices = new ArrayList<>();
        for(int i=0;i<t1.size();i++){
            if(((t1.getMeasurement(i, 0) - mean) / std) > z){
                indices.add(i);
            }
        }
        //Limit peaks based on taking only highest value of consecutive peaks
        ArrayList<Peak> finalPeaks = new ArrayList<>(indices.size());
        for(int i=0; i<indices.size(); i++){
            int max = indices.get(i);
            int furthest = i;
            for(int j=i+1;j<indices.size(); j++){
                if(indices.get(j-1)+1!=indices.get(j)) break;
                furthest = j;
                max = indices.get(j) > max ? indices.get(j) : max;
            }
            i = furthest;
            finalPeaks.add(new Peak(max, (float) t1.getMeasurement(max, 0)));
        }
        if(repPeak != null) finalPeaks = reducePeaks(finalPeaks, repPeak);
        return finalPeaks;
    }

    public static double[] calcFeatures(TimeSeries t1, Peak t1Peak) {
        return calcFeatures(t1, t1Peak, repTimeSeries);
    }

    /**
     *
     * @param bounds DetectedBounds representing a detected repetition for which we want to determine
     *               if it is a valid rep.
     * @return  True, if the result of logistic regression determines this is a valid rep, else false.
     */
    public static boolean accept(DetectedBounds bounds) {
        //TODO logistic regression to find coefficients
        final double b0 = 0, b1 = 1, b2 = 1, b3 = 1, b4 = 1, b5 = 1, b6 = 1;
        double res = b0 + (b1 * bounds.dst) + (b2 * bounds.max) + (b3 * bounds.min) + (b4 * bounds.sd)
                + (b5 * bounds.rms) + (b6 * bounds.dur);
        return 1 / (1 + Math.exp(-1.0 * res)) >= .5;
    }

    public static File loadRepetitionFile(String exerciseName, Context context) {
        File f = new File(context.getCacheDir(), exerciseName + "_calibration.dat");
        return f;
    }

    /**
     * Loads repPeak and repTimeSeries for this class from a reader containing three lines of data.
     * The first of which contains comma-separated values representing amplitudes at a consistent
     * time distance apart. The second line contains Peak information containing the index of the
     * peak in the stream, and the amplitude of the peak. The last line contains the index of the
     * major axis.
     * @param reader    reader from which to read the TimeSeries and Peak data.
     */
    public static void loadRepetitionPatternTimeSeries(BufferedReader reader) {
        TimeSeriesBase.Builder builder = TimeSeriesBase.builder();
        try {
            String line = reader.readLine();
            if (line == null)
                return;
            String[] numbers = line.split(",");
            int i = 0;
            for (String s : numbers) {
                try {
                    builder = builder.add(i++, Float.parseFloat(s));
                } catch (NumberFormatException e){
                    Log.d("DataUtils", "Format exception when reading from file");
                    return;
                }
            }

            line = reader.readLine();
            numbers = line.split(",");
            repPeak = new Peak(Integer.parseInt(numbers[0]), Float.parseFloat(numbers[1]));
            repTimeSeries = builder.build();

            line = reader.readLine();
            majorAxisIndex = Integer.parseInt(line);

            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Class representing a detected repetitions bounds and stats about the detected rep.
     * dst - DTW distance to the repetition pattern time series
     * max - maximal value of amplitude the stream takes on (also the Peak amplitude)
     * min - minimal value of amplitude the stream takes on
     * sd  - standard deviation of the repetition amplitude values
     * rms - root mean square of the repetition amplitude values
     * s   - start time of the repetition
     * e   - end time of the repetition
     * dur - duration of the repetition, calculated as e - s
     */
    static class DetectedBounds {
        double dst, max, min, sd, rms, dur;
        int s, e;

        public DetectedBounds(int s, int e, double dst, double max, double min, double sd, double rms) {
            this.s = s;
            this.e = e;
            this.dst = dst;
            this.max = max;
            this.min = min;
            this.sd = sd;
            this.rms = rms;
            this.dur = e - s;
        }
    }

    public static class Peak {
        int index;
        float amplitude;

        public Peak(int index, float amplitude) {
            this.index = index;
            this.amplitude = amplitude;
        }
    }

    public interface Listener {
        void respondToRep();
    }

}
