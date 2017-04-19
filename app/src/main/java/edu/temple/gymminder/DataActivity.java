package edu.temple.gymminder;

import android.content.Intent;
import android.net.Uri;
import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;

/**
 * Handles vibration on reps, calling DataUtils methods to handle data, accelerometer recording,
 * responding to calling component with relevant data such as number of reps, stream stats, the
 * stream array, etc.
 */
public class DataActivity extends AppCompatActivity implements
        ExerciseDataFragment.OnFragmentInteractionListener,
        CalibrateFragment.OnFragmentInteractionListener {

    public final static String EXTRA_REPS_DONE = "We smoked the last one an hour ago";
    public final static String EXTRA_MAX_VELOCITY = "Cathy I'm lost";
    public final static String EXTRA_AVG_VELOCITY = "I don't know why";

    private Exercise mExercise;
    private ArrayList<ArrayList<Float>> data;

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_data);

        if (savedInstanceState != null) {
            mExercise = (Exercise) savedInstanceState.getSerializable(DetailFragment.EXTRA_EXERCISE);
        } else {
            Bundle extras = getIntent().getExtras();
            if (extras == null)
                throw new RuntimeException("Attempting to create DataActivity without an exercise");

            mExercise = (Exercise) extras.getSerializable(DetailFragment.EXTRA_EXERCISE);
        }

        // Try to load the repetition pattern data for this exercise. Calibrate for it if it doesn't exist.
        File f = DataUtils.loadRepetitionFile(mExercise.name, this);
        if (f.exists()) {
            try {
                BufferedReader reader = new BufferedReader(new FileReader(f));
                DataUtils.loadRepetitionPatternTimeSeries(reader);
                goToExerciseData();
            } catch (FileNotFoundException e) {
                goToCalibrate();
            }

        } else {
            goToCalibrate();
        }
    }

    void result(int reps, int mv, int av) {
        Intent intent = new Intent();
        intent.putExtra(EXTRA_REPS_DONE, reps);
        intent.putExtra(EXTRA_MAX_VELOCITY, mv);
        intent.putExtra(EXTRA_AVG_VELOCITY, av);
        this.setResult(RESULT_OK, intent);
        finish();
    }

    void result(int reps) {
        Intent intent = new Intent();
        float[] f = DataUtils.maxAndAvg(DataUtils.riemann(data.get(2)));
        intent.putExtra(EXTRA_REPS_DONE, reps);
        intent.putExtra(EXTRA_MAX_VELOCITY, f[0]);
        intent.putExtra(EXTRA_AVG_VELOCITY, f[1]);
        this.setResult(RESULT_OK, intent);
        finish();
    }

    public void goToExerciseData() {
        getSupportFragmentManager().beginTransaction()
                .replace(R.id.data_frame, ExerciseDataFragment.newInstance(mExercise))
                .commit();
    }

    public void goToCalibrate() {
        getSupportFragmentManager().beginTransaction()
                .replace(R.id.data_frame, CalibrateFragment.newInstance(mExercise.name))
                .commit();
    }

    @Override
    public void onCalibrationComplete(Uri uri, Integer majorAxisIndex) {
        File f = new File(uri.getPath());
        if (f.exists()) {
            try {
                BufferedReader reader = new BufferedReader(new FileReader(f));
                DataUtils.loadRepetitionPatternTimeSeries(reader);
                goToExerciseData();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
    }

    @Override
    public void didFinish(int reps, int mv, int av, ArrayList<ArrayList<Float>> data) {
        result(reps, mv, av);
    }
}
