package edu.temple.gymminder;

import android.content.Context;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.net.Uri;
import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.view.LayoutInflater;
import android.view.Menu;
import android.view.MenuInflater;
import android.view.MenuItem;
import android.view.View;
import android.view.ViewGroup;
import android.widget.Button;
import android.widget.ProgressBar;
import android.widget.TextView;

import com.fastdtw.timeseries.TimeSeries;
import com.fastdtw.timeseries.TimeSeriesBase;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.temple.gymminder.processing.DataProcessor;
import edu.temple.gymminder.processing.DataUtils;

/**
 * A simple {@link Fragment} subclass.
 * Activities that contain this fragment must implement the
 * {@link CalibrateFragment.OnFragmentInteractionListener} interface
 * to handle interaction events.
 * Use the {@link CalibrateFragment#newInstance} factory method to
 * create an instance of this fragment.
 */
public class CalibrateFragment extends Fragment {

    private static final String EXERCISE_ARG = "ARG_EXERCISE";

    private enum ButtonState {
        START, STOP, REDO
    }

    private OnFragmentInteractionListener mListener;

    private SensorManager mSensorManager;
    private SensorEventListener mSensorListener;

    private String mExerciseName;
    private ArrayList<Float> xValues, yValues, zValues;
    private ArrayList<Long> timestamps;

    // UI
    private ProgressBar mProgressBar;
    private Button mButton;
    private ButtonState mButtonState = ButtonState.START;

    public CalibrateFragment() {
        // Required empty public constructor
    }

    /**
     * Use this factory method to create a new instance of
     * this fragment using the provided parameters.
     *
     * @param exerciseName Parameter 1.
     * @return A new instance of fragment CalibrateFragment.
     */
    // TODO: Rename and change types and number of parameters
    public static CalibrateFragment newInstance(String exerciseName) {
        CalibrateFragment fragment = new CalibrateFragment();
        Bundle args = new Bundle();
        args.putString(EXERCISE_ARG, exerciseName);
        fragment.setArguments(args);
        return fragment;
    }

    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);

        // Show "Done" button in menu when this fragment is present
        setHasOptionsMenu(true);

        if (getArguments() != null) {
            mExerciseName = getArguments().getString(EXERCISE_ARG);
        }
    }

    @Override
    public View onCreateView(LayoutInflater inflater, ViewGroup container,
                             Bundle savedInstanceState) {
        // Inflate the layout for this fragment
        View v = inflater.inflate(R.layout.fragment_calibrate, container, false);

        xValues = new ArrayList<>();
        yValues = new ArrayList<>();
        zValues = new ArrayList<>();
        timestamps = new ArrayList<>();

        TextView titleTextView = (TextView) v.findViewById(R.id.calibrate_title);
        titleTextView.setText(mExerciseName);

        mProgressBar = (ProgressBar) v.findViewById(R.id.calibrate_progress);
        mButton = (Button) v.findViewById(R.id.calibrate_button);
        if (mButton != null) {
            mButton.setOnClickListener(new View.OnClickListener() {
                @Override
                public void onClick(View view) {
                    onButtonPressed();
                }
            });
        }

        return v;
    }

    @Override
    public void onCreateOptionsMenu(Menu menu, MenuInflater inflater) {
        inflater.inflate(R.menu.fragment_calibrate_menu, menu);
    }

    public void onButtonPressed() {

        switch (mButtonState) {
            case REDO:
                xValues.clear();
                yValues.clear();
                zValues.clear();
                timestamps.clear();
            case START:
                mButton.setText(R.string.done_calibration);
                mButtonState = ButtonState.STOP;

                // Animate and show the progress bar
                mProgressBar.animate();
                mProgressBar.setVisibility(View.VISIBLE);

                setupSensor();
                break;
            case STOP:
                mButton.setText(R.string.redo_calibration);
                mButtonState = ButtonState.REDO;

                // Hide the progress bar
                mProgressBar.setVisibility(View.INVISIBLE);

                mSensorManager.unregisterListener(mSensorListener);
                break;
            default:
                break;
        }
    }

    void setupSensor() {
        mSensorManager = (SensorManager) getContext().getSystemService(Context.SENSOR_SERVICE);
        Sensor sensor = mSensorManager.getDefaultSensor(Sensor.TYPE_LINEAR_ACCELERATION);

        mSensorListener = new SensorEventListener() {
            @Override
            public void onSensorChanged(SensorEvent event) {
                xValues.add(event.values[0]);
                yValues.add(event.values[1]);
                zValues.add(event.values[2]);
                timestamps.add(event.timestamp);
            }

            @Override
            public void onAccuracyChanged(Sensor sensor, int accuracy) {}
        };

        mSensorManager.registerListener(mSensorListener, sensor, (int) DataProcessor.POLLING_RATE);
    }

    void process() {
        File f = DataUtils.loadRepetitionFile(mExerciseName, getContext());
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(f));
            StringBuilder sb = new StringBuilder();
            ArrayList<ArrayList<Float>> axes = new ArrayList<>(3);
            axes.add(xValues);
            axes.add(yValues);
            axes.add(zValues);

            // Find the major axis among the three
            int majorAxisIndex = DataUtils.detectMajorAxis(axes);
            List<Float> filtered = DataUtils.applySavitzkyGolayFilter(axes.get(majorAxisIndex));

            // Build a time series to use for peak detection
            TimeSeriesBase.Builder builder = TimeSeriesBase.builder();
            int i = 0;
            for (Float val : filtered) {
                builder.add(i++, val);
            }
            TimeSeries timeSeries = builder.build();

            // Find the peak of the rep
            ArrayList<DataUtils.Peak> peaks = DataUtils.zScorePeakDetection(timeSeries);
            DataUtils.Peak peak;
            if (peaks.size() > 0) {
                // use the peak with the max value
                peak = peaks.get(0);
                for (DataUtils.Peak p : peaks) {
                    if (p.amplitude > peak.amplitude)
                        peak = p;
                }
            } else {
                // TODO: Tell the user to redo the rep, it wasn't good enough to find a peak
                return;
            }

            // Write first line (amplitudes)
            for (Float val : filtered) {
                sb.append(val);
                sb.append(",");
            }
            writer.append(sb.toString());
            writer.newLine();

            // Write second line (peak info)
            writer.write(peak.index + "," + peak.amplitude);
            writer.newLine();

            // Write third line (index of major axis)
            writer.write(majorAxisIndex + "");

            writer.close();
            if (mListener != null) {
                mListener.onCalibrationComplete(Uri.fromFile(f), majorAxisIndex);
            }
            return;
        } catch (IOException e) {
            e.printStackTrace();
        }

        mListener.onCalibrationComplete(null, -1);
    }

    @Override
    public void onAttach(Context context) {
        super.onAttach(context);
        if (context instanceof OnFragmentInteractionListener) {
            mListener = (OnFragmentInteractionListener) context;
        } else {
            throw new RuntimeException(context.toString()
                    + " must implement OnFragmentInteractionListener");
        }
    }

    @Override
    public void onDetach() {
        super.onDetach();
        mListener = null;
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        // handle item selection
        switch (item.getItemId()) {
            case R.id.done_option:
                if (mButtonState == ButtonState.REDO)
                    process();
                return true;
            default:
                return super.onOptionsItemSelected(item);
        }
    }

    public interface OnFragmentInteractionListener {

        /**
         * Called when the fragment has finished the calibration process.
         * @param uri The uri to the file containing the data of the processed time series recorded during calibration
         * @param majorAxisIndex The index of the detected major axis
         */
        void onCalibrationComplete(Uri uri, Integer majorAxisIndex);
    }
}
