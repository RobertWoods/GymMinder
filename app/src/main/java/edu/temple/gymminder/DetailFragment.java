package edu.temple.gymminder;


import android.app.Activity;
import android.content.Context;
import android.content.Intent;
import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.util.Log;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.BaseAdapter;
import android.widget.Button;
import android.widget.LinearLayout;
import android.widget.ListView;
import android.widget.TextView;

import com.google.firebase.auth.FirebaseAuth;

import java.util.ArrayList;
import java.util.Calendar;

import edu.temple.gymminder.models.Exercise;
import edu.temple.gymminder.models.Workout;


/**
 * A Fragment showing a stored workout. This contains information about the workout, including
 * the list of exercises defined in the workout, the current number of sets out of the total number
 * of sets the user has completed, and a button that takes the user to DataActivity to allow them
 * to record their set.
 */
public class DetailFragment extends Fragment {

    private static final String WORKOUT_ARG = "SPARROW THAN A SNAIL";
    private static final String NAME_ARG = "There are no registered users";
    private static final String PLACEHOLDER_WORKOUT = "So for now let's just";
    public static final String EXTRA_EXERCISE = "I heard cathedral bells";
    public static final int RESULT_REPS = 7073;

    ListView lv;
    String workoutName;

    public DetailFragment() {
        // Required empty public constructor
    }

    public static DetailFragment newInstance(Workout workout, String name) {
        DetailFragment detailFragment = new DetailFragment();
        Bundle args = new Bundle();
        args.putSerializable(WORKOUT_ARG, workout);
        detailFragment.setArguments(args);
        args.putString(NAME_ARG, name);
        return detailFragment;
    }


    @Override
    public View onCreateView(LayoutInflater inflater, ViewGroup container,
                             Bundle savedInstanceState) {
        // Inflate the layout for this fragment
        View v = inflater.inflate(R.layout.fragment_detail, container, false);
        Workout workout = (Workout) getArguments().getSerializable(WORKOUT_ARG);
        workoutName = getArguments().getString(NAME_ARG);
        if (workout == null) return v;
        ((TextView) v.findViewById(R.id.textView)).setText(
                workoutName != null ? workoutName : workout.toString());

        final ArrayList<Exercise> exercises = workout.exercises;
        lv = (ListView) v.findViewById(R.id.workoutsList);
        lv.setAdapter(new BaseAdapter() {
            @Override
            public int getCount() {
                return exercises.size();
            }

            @Override
            public Object getItem(int position) {
                return exercises.get(position);
            }

            @Override
            public long getItemId(int position) {
                return position;
            }

            @Override
            public View getView(final int position, View convertView, ViewGroup parent) {
                //TODO refactor into using Fragments maybe
                LayoutInflater inflater = (LayoutInflater) getContext().getSystemService(Context.LAYOUT_INFLATER_SERVICE);
                final Exercise exercise = (Exercise) getItem(position);
                if (exercise.completed == null) exercise.initActive();
                View item = inflater.inflate(R.layout.item_exercise, parent, false); //TODO view reuse
                ((TextView) item.findViewById(R.id.workoutName)).setText(exercise.toString());
                final TextView setProgress = (TextView) item.findViewById(R.id.setNumber);
                setProgress.setText(getString(R.string.sets_progress, 0, exercise.sets));
                Button startButton = (Button) item.findViewById(R.id.startTrackerButton);
                if (exercise.setsDone >= exercise.sets)
                    startButton.setVisibility(View.INVISIBLE);
                else {
                    startButton.setOnClickListener(new View.OnClickListener() {
                        @Override
                        public void onClick(View v) {
                            setProgress.setText(getString(R.string.sets_progress,
                                    ++exercise.setsDone, exercise.sets));
                            Intent intent = new Intent(getContext(), DataActivity.class);
                            intent.putExtra(EXTRA_EXERCISE, exercise);
                            startActivityForResult(intent, RESULT_REPS);
                        }
                    });
                }
                return item;
            }
        });
        return v;
    }

    @Override
    public void onActivityResult(int requestCode, int resultCode, Intent data) {
        if (requestCode == RESULT_REPS) {
            if (resultCode == Activity.RESULT_OK) {
                int repsDone = data.getIntExtra(DataActivity.EXTRA_REPS_DONE, -1);
                Log.d("OuO", String.valueOf(repsDone));
                for (int i = 0; i < lv.getAdapter().getCount(); i++) {
                    Exercise ex = (Exercise) lv.getAdapter().getItem(i);
                    if (ex.setsDone != ex.sets) return;
                }
                onWorkoutFinished();
            }
        } else {
            Log.d("OnO", "Something went wrong");
        }
    }

    public void onWorkoutFinished() {
        //Add button to let user save workout to database, or maybe just do it automatically
        Button button = new Button(getContext());
        button.setText(R.string.save_workout_text);
        button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                DbHelper dbHelper = DbHelper.newInstance(null);
                ArrayList<Exercise> exercises = new ArrayList<>();
                for (int i = 0; i < lv.getAdapter().getCount(); i++) {
                    exercises.add((Exercise) lv.getAdapter().getItem(i));
                }
                dbHelper.addWorkout(new Workout(exercises),
                        workoutName != null ? workoutName : PLACEHOLDER_WORKOUT,
                        FirebaseAuth.getInstance().getCurrentUser(),
                        Calendar.getInstance().getTime());
            }
        });
        ((LinearLayout) getView()).addView(button);
    }

}
