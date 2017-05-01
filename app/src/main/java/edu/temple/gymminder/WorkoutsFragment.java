package edu.temple.gymminder;


import android.content.Context;
import android.os.Bundle;
import android.support.v4.app.Fragment;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.ListView;
import android.widget.TabHost;
import android.widget.TextView;

import com.google.firebase.auth.FirebaseAuth;

import java.util.ArrayList;
import java.util.Map;

import edu.temple.gymminder.models.Exercise;
import edu.temple.gymminder.models.Workout;


/**
 * Currently handles workout selection, and navigation to workout creator. This is the Fragment shown
 * after the user has been authenticated.
 */
public class WorkoutsFragment extends Fragment implements DbHelper.Listener {

    FirebaseAuth auth = FirebaseAuth.getInstance();
    DbHelper db = DbHelper.newInstance(this);
    DetailListener listener;


    public WorkoutsFragment() {
        // Required empty public constructor
    }


    @Override
    public View onCreateView(LayoutInflater inflater, ViewGroup container,
                             Bundle savedInstanceState) {
        // Inflate the layout for this fragment
        View v = inflater.inflate(R.layout.fragment_workouts, container, false);

        TabHost host = (TabHost) v.findViewById(R.id.tabHost);
        host.setup();

        //Tab 1
        TabHost.TabSpec spec = host.newTabSpec(getString(R.string.title_workouts));
        spec.setContent(R.id.workoutsTab);
        spec.setIndicator(getString(R.string.title_workouts));
        host.addTab(spec);

        //Tab 2
        spec = host.newTabSpec(getString(R.string.title_exercises));
        spec.setContent(R.id.exercisesTab);
        spec.setIndicator(getString(R.string.title_exercises));
        host.addTab(spec);

        v.findViewById(R.id.add_workout_fab).setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                listener.goToWorkoutCreator();
            }
        });

        ((TextView) v.findViewById(R.id.greeting2)).setText(
                getResources().getString(R.string.greeting2 ));
        return v;
    }

    @Override
    public void onStart() {
        super.onStart();
        db.retrieveAllWorkouts(auth.getCurrentUser());
        db.getCatalog();
    }

    @Override
    public void onResume(){
        super.onResume();
    }

    @Override
    public void updateUi(Workout workout) {}

    @Override
    public void respondToWorkouts(final ArrayList<Workout> workouts, final ArrayList<String> names) {
        ((TextView) getView().findViewById(R.id.workouts)).setText(
                getResources().getText(R.string.stored_workouts_text));
        ListView lv = (ListView) getView().findViewById(R.id.workoutsList);

        ArrayAdapter<String> adapter = new ArrayAdapter<>(getActivity(), android.R.layout.simple_list_item_1, names);
        lv.setAdapter(adapter);
        lv.setOnItemClickListener(new AdapterView.OnItemClickListener() {
            @Override
            public void onItemClick(AdapterView<?> parent, View view, int position, long id) {
                listener.goToWorkoutsDetail(workouts.get(position), (String) parent.getAdapter().getItem(position));
            }
        });
    }

    @Override
    public void respondToHistory(ArrayList<Workout> workouts, ArrayList<String> names, ArrayList<String> workoutNames, Map<String, String> dates) {

    }

    @Override
    public void respondToCatalog(final ArrayList<Exercise> exercises) {
        ListView lv = (ListView) getView().findViewById(R.id.exerciseListView);

        final ArrayList<String> names = new ArrayList<>(exercises.size());
        for (Exercise e : exercises) {
            names.add(e.name);
        }
        ArrayAdapter<String> adapter = new ArrayAdapter<>(getActivity(), android.R.layout.simple_list_item_1, names);
        lv.setAdapter(adapter);
        lv.setOnItemClickListener(new AdapterView.OnItemClickListener() {
            @Override
            public void onItemClick(AdapterView<?> parent, View view, int position, long id) {
                if (exercises.get(position) != null)
                    listener.goToAdHocCreator(exercises, position);
            }
        });
    }

    @Override
    public void onWorkoutAdded() {
        ListView lv = (ListView) getView().findViewById(R.id.workoutsList);
        ArrayAdapter aa = (ArrayAdapter) lv.getAdapter();
        aa.notifyDataSetChanged();
    }

    @Override
    public void onAttach(Context c) {
        super.onAttach(c);
        listener = (DetailListener) c;
    }

    public interface DetailListener {
        void goToWorkoutsDetail(Workout workout, String name);
        void goToWorkoutCreator();
        void goToAdHocCreator(ArrayList<Exercise> exercises, int selected);
    }

}
