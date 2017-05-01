package edu.temple.gymminder;

import android.content.Context;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.os.Bundle;
import android.os.Vibrator;
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

import java.util.ArrayList;

import edu.temple.gymminder.models.Exercise;

import static android.content.Context.VIBRATOR_SERVICE;

/**
 * A simple {@link Fragment} subclass.
 * Activities that contain this fragment must implement the
 * {@link ExerciseDataFragment.OnFragmentInteractionListener} interface
 * to handle interaction events.
 * Use the {@link ExerciseDataFragment#newInstance} factory method to
 * create an instance of this fragment.
 */
public class ExerciseDataFragment extends Fragment implements DataUtils.Listener {

    private static final String ARG_EXERCISE = "ARG_EXERCISE";

    private Exercise mExercise;
    private ArrayList<ArrayList<Float>> data = new ArrayList<>(4);
    private ArrayList<Long> timestamps = new ArrayList<>();
    private int mReps = 0;
    private Vibrator vibrator;

    private SensorManager mSensorManager;
    private SensorEventListener mSensorListener;

    private TextView mRepsTextView;
    private ProgressBar mSetProgressBar;

    private OnFragmentInteractionListener mListener;

    public ExerciseDataFragment() {
        // Required empty public constructor
    }

    /**
     * Use this factory method to create a new instance of
     * this fragment using the provided parameters.
     *
     * @param exercise The current exercise
     * @return A new instance of fragment ExerciseDataFragment.
     */
    // TODO: Rename and change types and number of parameters
    public static ExerciseDataFragment newInstance(Exercise exercise) {
        ExerciseDataFragment fragment = new ExerciseDataFragment();
        Bundle args = new Bundle();
        args.putSerializable(ARG_EXERCISE, exercise);
        fragment.setArguments(args);
        return fragment;
    }

    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        if (getArguments() != null) {
            mExercise = (Exercise) getArguments().getSerializable(ARG_EXERCISE);
        }

        // Show "Done" button in menu when this fragment is present
        setHasOptionsMenu(true);

        vibrator = (Vibrator) getActivity().getSystemService(VIBRATOR_SERVICE);
        for (int i = 0; i < 3; i++) data.add(new ArrayList<Float>());
    }

    @Override
    public View onCreateView(LayoutInflater inflater, ViewGroup container,
                             Bundle savedInstanceState) {
        // Inflate the layout for this fragment
        View v = inflater.inflate(R.layout.fragment_exercise_data, container, false);

        Button button = (Button) v.findViewById(R.id.finish_set_button);
        button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                onButtonPressed();
            }
        });

        ((TextView) v.findViewById(R.id.exercise_title)).setText(mExercise.name);

        mRepsTextView = (TextView) v.findViewById(R.id.rep_text_view);
        mSetProgressBar = (ProgressBar) v.findViewById(R.id.setProgressBar);

        return v;
    }

    @Override
    public void onCreateOptionsMenu(Menu menu, MenuInflater inflater) {
        inflater.inflate(R.menu.fragment_exercise_data_menu, menu);
    }

    public void onButtonPressed() {
        if (mListener != null) {
            mListener.didFinish(mReps, data.get(DataUtils.majorAxisIndex));
        }
    }

    void setupSensor() {
        mSensorManager = (SensorManager) getActivity().getSystemService(Context.SENSOR_SERVICE);
        mSensorListener = new SensorEventListener() {
            @Override
            public void onSensorChanged(SensorEvent event) {
                DataUtils.process(event.values, event.timestamp);
            }

            @Override
            public void onAccuracyChanged(Sensor sensor, int accuracy) {}
        };
        Sensor sensor = mSensorManager.getDefaultSensor(Sensor.TYPE_LINEAR_ACCELERATION);
        DataUtils.init(data, timestamps);
        DataUtils.setListener(this);
        mSensorManager.registerListener(mSensorListener, sensor, (int) DataUtils.POLLING_RATE);
    }

    private void finish() {
        mSensorManager.unregisterListener(mSensorListener);
        DataUtils.removeListener();
        mListener.didFinish(mReps, data.get(DataUtils.majorAxisIndex));
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
    public void onResume() {
        super.onResume();
        setupSensor();
    }

    @Override
    public void onPause() {
        super.onPause();
        DataUtils.removeListener();
        mSensorManager.unregisterListener(mSensorListener);
    }

    @Override
    public void respondToRep() {

        getActivity().runOnUiThread(new Runnable() {
            @Override
            public void run() {
                if (vibrator.hasVibrator())
                    vibrator.vibrate(100);

                mReps++;
                mRepsTextView.setText(mReps + "");

                if (mExercise != null) {
                    int progress = (int) ((mReps / (double) mExercise.reps) * 100);
                    mSetProgressBar.setProgress(progress);

                    // We're done!
                    if (mExercise.reps == mReps)
                        finish();
                }
            }
        });

    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        // handle item selection
        switch (item.getItemId()) {
            case R.id.recalibrate_option:
                mListener.requestRecalibration();
                return true;
            default:
                return super.onOptionsItemSelected(item);
        }
    }

    public interface OnFragmentInteractionListener {
        void didFinish(int reps, ArrayList<Float> data);
        void requestRecalibration();
    }
}
