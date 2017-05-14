package edu.temple.gymminder.processing;

import android.util.Log;
import android.util.SparseArray;

import com.fastdtw.timeseries.TimeSeries;
import com.fastdtw.timeseries.TimeSeriesBase;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import edu.temple.gymminder.models.Peak;

/**
 * Created by rober_000 on 5/13/2017.
 */

public class DataProcessor {

    public static final float MS2S_CONVERSION = 1.0f / 1000000.0f;
    private static final long ERROR = 1000;
    public static final long SECOND = 1000000;
    public static final long POLLING_FREQUENCY = 30;
    public static final long POLLING_RATE = SECOND / POLLING_FREQUENCY;

    private float[] avgNode = null;
    private ArrayList<ArrayList<Float>> data;
    private ArrayList<ArrayList<Float>> processedData;
    private ArrayList<Long> timestamps;

    public SparseArray<Peak> peaks = new SparseArray<>();
    public Peak repPeak;
    public TimeSeries repTimeSeries;
    public int majorAxisIndex;

    private ExecutorService executorService;
    private DataUtils.Listener listener;

    DataProcessor(ArrayList<ArrayList<Float>> dataList, ArrayList<Long> time) {
        data = dataList;
        timestamps = time;
        processedData = new ArrayList<>(3);
        peaks = new SparseArray<>();
        for (int i = 0; i < 3; i++) processedData.add(new ArrayList<Float>());
    }

    DataProcessor(ArrayList<ArrayList<Float>> dataList, ArrayList<Long> time, ArrayList<ArrayList<Float>> processedDataList) {
        this(dataList, time);
        processedData = processedDataList;
        avgNode = null;
    }

    DataProcessor(ArrayList<ArrayList<Float>> dataList, ArrayList<Long> time, ArrayList<ArrayList<Float>> processedDataList, File file) {
        this(dataList, time, processedDataList);
        try {
            loadRepetitionPatternTimeSeries(new BufferedReader(new FileReader(file)));
            peaks = new SparseArray<>((int) (DataUtils.EXPANSION_VALUE * (repTimeSeries.size() - repPeak.index)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }

    void setListener(DataUtils.Listener l) {
        listener = l;
        executorService = Executors.newSingleThreadExecutor();
    }

    //Need this to prevent possible memory leak
    void removeListener() {
        executorService.shutdownNow();
        listener = null;
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
     void addToProcessQueue(final float[] values, final long timestamp) {
        executorService.submit(new Runnable() {
            @Override
            public void run() {
                process(values, timestamp);
            }
        });

    }

    void process(float[] values, long timestamp) {
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
        long longDuration = timestamp - timestamps.get(timestamps.size() - 1);

        if ((longDuration + ERROR) < POLLING_RATE || avgNode != null) {
            //average the points with sum node
            avgNode = DataUtils.average(avgNode, x, duration);
            duration = avgNode[1];
            x = avgNode[0];
        }
        if ((longDuration + ERROR) >= (1.0 / POLLING_FREQUENCY)) {
            //We can approximate timestamp value by adding .1s to previous value
            //Maybe not the best idea since it (maybe) causes drift when we interpolate, idk :d
            timestamps.add(timestamps.get(timestamps.size() - 1) + POLLING_RATE);
            data.get(i).add(x);
            int size = processedData.get(i).size();
            for (int j = size; j >= size - (DataUtils.SG_FILTER.length / 2) && j > 0; j--) {
                /*
                    We want to re-process any data points that didn't have enough data to the right
                    for the entire filter to run on. The alternative to this is to delay the signal
                    by waiting until we have enough data points for the entire filter.
                 */
                DataUtils.applySGFilterRealtime(j, data.get(i), processedData.get(i));
            }
            avgNode = null;
            Peak newPeak = DataUtils.detectPeak(size, processedData.get(majorAxisIndex), repPeak, 1);
            //TODO: Probably want to put this in a thread that enqueues new peaks to check
            if (newPeak != null) {
                    /*
                        The key for our HashMap is the index of the new peak, which is the current
                        size of processedData at the time of insertion, plus the difference between
                        the size of repTimeSeries and the index of repPeak multiplied by
                        EXPANSION_VALUE.
                     */
                peaks.put((int) (newPeak.index + DataUtils.EXPANSION_VALUE *
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
                DataUtils.DetectedBounds bounds = DataUtils.detectBounds(t1, peaks.get(processedData.get(i).size()), repTimeSeries, repPeak);
                int peakIndex = peaks.get(processedData.get(i).size()).index;
                peaks.remove(processedData.get(i).size());
                if (DataUtils.accept(bounds)) {
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

    /**
     * Loads repPeak and repTimeSeries for this class from a reader containing three lines of data.
     * The first of which contains comma-separated values representing amplitudes at a consistent
     * time distance apart. The second line contains Peak information containing the index of the
     * peak in the stream, and the amplitude of the peak. The last line contains the index of the
     * major axis.
     *
     * @param reader reader from which to read the TimeSeries and Peak data.
     */
    public void loadRepetitionPatternTimeSeries(BufferedReader reader) {
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
                } catch (NumberFormatException e) {
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
}
