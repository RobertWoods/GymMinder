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

import com.google.firebase.auth.FirebaseAuth;

import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.Map;

import edu.temple.gymminder.models.Exercise;
import edu.temple.gymminder.models.Workout;


/**
 * A simple {@link Fragment} subclass.
 * Activities that contain this fragment must implement the
 * {@link HistoryFragment.OnFragmentInteractionListener} interface
 * to handle interaction events.
 */
public class HistoryFragment extends Fragment implements DbHelper.Listener {




    public ArrayList<String> workoutsNames2;
    ArrayList<Workout> workouts2;
    public int q = 0;
    public boolean q1 = false;
    FirebaseAuth auth = FirebaseAuth.getInstance();
    DbHelper db = DbHelper.newInstance(this);
    WorkoutsFragment.DetailListener listener;
    private OnFragmentInteractionListener mListener;

    public HistoryFragment() {
        // Required empty public constructor
    }

    @Override
    public View onCreateView(LayoutInflater inflater, ViewGroup container,
                             Bundle savedInstanceState) {
        // Inflate the layout for this fragment
        View view = inflater.inflate(R.layout.fragment_history, container, false);
        return view;
    }

    @Override
    public void onAttach(Context context) {
        super.onAttach(context);
        if (context instanceof OnFragmentInteractionListener) {
            mListener = (OnFragmentInteractionListener) context;
            listener = (WorkoutsFragment.DetailListener) context;
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
    public void updateUi(Workout workout) {

    }

    @Override
    public void respondToWorkouts(ArrayList<Workout> workouts, ArrayList<String> names) {

    }

    @Override
    public void onResume(){
        super.onResume();
        db.getHistory(auth.getCurrentUser());
    }

    @Override
    public void respondToHistory(final ArrayList<Workout> workouts, final ArrayList<String> names, final ArrayList<String> workoutNames, final Map<String, String> dates) {
        ListView lv = (ListView) getView().findViewById(R.id.historylist);
        ArrayList<String> names2 = new ArrayList<String>();
        for(String x : names){
            String[] y = x.split("\\s+");
            int day = Integer.parseInt(y[0]);
            int year = Integer.parseInt(y[1]);
            final Calendar calendar = Calendar.getInstance();
            calendar.set(Calendar.DAY_OF_YEAR, day);
            final Date date = new Date(calendar.getTimeInMillis());
            Calendar cal = Calendar.getInstance();
            cal.setTime(date);
            int year2 = cal.get(Calendar.YEAR);
            int month = cal.get(Calendar.MONTH);
            int day2 = cal.get(Calendar.DAY_OF_MONTH);
            month++;
            names2.add(month + "/" + day2 + "/" + year);
        }


        q = -1;
        q1 = false;
        ArrayAdapter<String> adapter = new ArrayAdapter<>(getActivity(), android.R.layout.simple_list_item_1, names2);

        lv.setAdapter(adapter);
        lv.setOnItemClickListener(new AdapterView.OnItemClickListener() {
            @Override
            public void onItemClick(AdapterView<?> parent, View view, int position, long id) {
                ListView lv = (ListView) getView().findViewById(R.id.historylist);
                if(q1){
                    listener.goToWorkoutsDetail(workouts2.get(position), workoutsNames2.get(position));
                }
                else{


                workoutsNames2 = new ArrayList<String>();
                workouts2 = new ArrayList<Workout>();
                //for (Workout workout1 : workouts){
                    for(int i = 0; i<workouts.size(); i++){
                        if(names.get(position).equals(dates.get(workoutNames.get(i)))){
                            workoutsNames2.add(workoutNames.get(i));
                            workouts2.add(workouts.get(i));
                        }

                }



            }
            ArrayAdapter<String> adapter = new ArrayAdapter<>(getActivity(), android.R.layout.simple_list_item_1, workoutsNames2);
            q1 = true;
            lv.setAdapter(adapter);
            }
        });
    }

    @Override
    public void respondToCatalog(ArrayList<Exercise> exercises) {}

    @Override
    public void onWorkoutAdded() {

    }

    /**
     * This interface must be implemented by activities that contain this
     * fragment to allow an interaction in this fragment to be communicated
     * to the activity and potentially other fragments contained in that
     * activity.
     * <p>
     * See the Android Training lesson <a href=
     * "http://developer.android.com/training/basics/fragments/communicating.html"
     * >Communicating with Other Fragments</a> for more information.
     */
    public interface OnFragmentInteractionListener {
    }
}
