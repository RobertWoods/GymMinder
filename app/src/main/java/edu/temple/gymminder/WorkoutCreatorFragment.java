package edu.temple.gymminder;


import android.content.Context;
import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.text.Editable;
import android.text.TextWatcher;
import android.util.TypedValue;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.AdapterView;
import android.widget.BaseAdapter;
import android.widget.EditText;
import android.widget.LinearLayout;
import android.widget.ListView;
import android.widget.Spinner;
import android.widget.TextView;

import com.google.firebase.auth.FirebaseAuth;

import java.util.ArrayList;

import edu.temple.gymminder.models.Exercise;
import edu.temple.gymminder.models.Workout;


/**
 * A simple {@link Fragment} subclass.
 */
public class WorkoutCreatorFragment extends Fragment {

    BaseAdapter listAdapter;
    ArrayList<Exercise> exercises = new ArrayList<>();
    String[] exerciseNames;
    Listener listener;
    String PLACEHOLDER_STRING = "You may be surprised to hear that when the revolution happens," +
            "the proletariat will be destroyed. Indeed, just as the bourgeois will be stripped" +
            "of power and erased, so too will the proletariat. After the revolution, there will" +
            "be no proletariat, and there will be no bourgeois. There will only be people." +
            "So rise up, comrades! And take what is yours!";
    //TODO maybe refactor this to something that makes more sense o3o
    WorkoutCreatorFragment self = this;

    public WorkoutCreatorFragment() {
        // Required empty public constructor
    }


    @Override
    public View onCreateView(LayoutInflater inflater, ViewGroup container,
                             Bundle savedInstanceState) {
        // Inflate the layout for this fragment
        View v = inflater.inflate(R.layout.fragment_workout_creator, container, false);
        listAdapter = new BaseAdapter() {
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
                //TODO view reuse
                LayoutInflater inflater = (LayoutInflater) getContext().getSystemService(
                        Context.LAYOUT_INFLATER_SERVICE);
                LinearLayout ll = (LinearLayout) inflater.inflate(R.layout.item_exercise_creator,
                        parent, false);

                final Spinner spinner = (Spinner) ll.findViewById(R.id.exerciseSpinner);
                spinner.setAdapter(new BaseAdapter() {
                    @Override
                    public int getCount() {
                        return exerciseNames.length;
                    }

                    @Override
                    public Object getItem(int i) {
                        return exerciseNames[i];
                    }

                    @Override
                    public long getItemId(int i) {
                        return i;
                    }

                    @Override
                    public View getView(int i, View view, ViewGroup viewGroup) {
                        LinearLayout ll = new LinearLayout(getContext());
                        TextView textView = new TextView(getContext());
                        textView.setText(getItem(i).toString());
                        //TODO make this level with EditTexts
                        textView.setTextSize(TypedValue.COMPLEX_UNIT_DIP, 18);
                        ll.addView(textView);
                        return ll;
                    }
                });
                spinner.setOnItemSelectedListener(new AdapterView.OnItemSelectedListener() {
                    @Override
                    public void onItemSelected(AdapterView<?> adapterView, View view, int i, long l) {
                        Exercise ex = (Exercise) listAdapter.getItem(position); //refers to ListView
                        String name = (String) spinner.getAdapter().getItem(i); //refers to Spinner
                        ex.setName(name);
                    }

                    @Override
                    public void onNothingSelected(AdapterView<?> adapterView) {

                    }
                });

                final Exercise exercise = (Exercise) getItem(position);
                if (!exercise.name.equals(PLACEHOLDER_STRING)) {
                    for (int i = 0; i < exerciseNames.length; i++) {
                        if (exercise.name.equals(exerciseNames[i])) {
                            spinner.setSelection(i);
                            break;
                        }
                    }
                }
                EditText editTextSets = (EditText) ll.findViewById(R.id.setsEditText);
                editTextSets.addTextChangedListener(new TextWatcher() {
                    @Override
                    public void beforeTextChanged(CharSequence charSequence, int i, int i1, int i2) {

                    }

                    @Override
                    public void onTextChanged(CharSequence charSequence, int i, int i1, int i2) {

                    }

                    @Override
                    public void afterTextChanged(Editable editable) {
                        try {
                            exercise.sets = Integer.parseInt(editable.toString());
                        } catch (Exception e) {
                            //TODO error checking
                        }
                    }
                });
                EditText editTextReps = (EditText) ll.findViewById(R.id.repsEditText);
                editTextReps.addTextChangedListener(new TextWatcher() {
                    @Override
                    public void beforeTextChanged(CharSequence charSequence, int i, int i1, int i2) {

                    }

                    @Override
                    public void onTextChanged(CharSequence charSequence, int i, int i1, int i2) {

                    }

                    @Override
                    public void afterTextChanged(Editable editable) {
                        try {
                            exercise.reps = Integer.parseInt(editable.toString());
                        } catch (Exception e) {
                            //TODO error checking
                        }
                    }
                });
                if (exercise.sets != -999) {
                    editTextSets.setText(String.valueOf(exercise.sets));
                }
                if (exercise.reps != -999) {
                    ((EditText) ll.findViewById(R.id.repsEditText)).setText(String.valueOf(exercise.reps));
                }
                return ll;
            }
        };

        ((ListView) v.findViewById(R.id.exercisesListview)).setAdapter(listAdapter);
        v.findViewById(R.id.addExerciseButton).setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                exercises.add(new Exercise(PLACEHOLDER_STRING, -999, -999));
                listAdapter.notifyDataSetChanged();
            }
        });
        v.findViewById(R.id.finishExerciseButton).setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                DbHelper dbHelper = DbHelper.newInstance(null);
                //TODO error checking for no input values in one of the exercise fields
                dbHelper.addNewWorkout(new Workout(exercises, true),
                        ((EditText) getView().findViewById(R.id.workoutNameEditText)).getText().toString(),
                        FirebaseAuth.getInstance().getCurrentUser());
                listener.finishFragment(self);
            }
        });
        return v;
    }

    @Override
    public void onAttach(Context c) {
        super.onAttach(c);
        exerciseNames = getResources().getStringArray(R.array.supported_exercises);
        listener = (Listener) c;
    }

    public interface Listener {
        void finishFragment(Fragment f);
    }

}
