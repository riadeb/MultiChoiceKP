import java.lang.reflect.Array;
import java.util.*;
import  java.io.File;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;
import  scpsolver.constraints.*;

import javax.rmi.ssl.SslRMIClientSocketFactory;


public class Reader {
    int p;
    int M;
    int K;
    int N;
    String filename;
    ArrayList<pair>[] data;
    ArrayList<pair>[] LP_filtered_data;
    int[] maxes_r;
    int upper_bnd_Rate;

    public Reader(String filename) throws Exception {
        File file = new File(filename);
        Scanner sc = new Scanner(file);
        this.N = Double.valueOf(sc.nextLine()).intValue();
        this.M = Double.valueOf(sc.nextLine()).intValue();
        this.K = Double.valueOf(sc.nextLine()).intValue();
        this.p = Double.valueOf(sc.nextLine()).intValue();
        this.filename = filename;
        data = new ArrayList[N];
        maxes_r = new int[N];
        LP_filtered_data = new ArrayList[N];


        int[][][] power_values = new int[N][K][M];
        int[][][] r_values = new int[N][K][M];

        for (int i = 0; i < N; i++) {
            data[i] = new ArrayList<pair>();
            for (int k = 0; k < K; k++) {
                String[] newline = sc.nextLine().stripLeading().split("   ");
                for (int m = 0; m < M; m++) {

                    power_values[i][k][m] = Double.valueOf(newline[m]).intValue();

                }
            }
        }
        for (int i = 0; i < N; i++) {
            int maxr = 0;
            for (int k = 0; k < K; k++) {
                String[] newline = sc.nextLine().stripLeading().split("   ");
                for (int m = 0; m < M; m++) {
                    r_values[i][k][m] = Double.valueOf(newline[m]).intValue();
                    maxr = Math.max(maxr,Double.valueOf(newline[m]).intValue());
                }
            }
            for(int j = 0; j <= i; j++) maxes_r[j] += maxr;
            this.upper_bnd_Rate += maxr;
        }
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < K; k++) {
                for (int m = 0; m < M; m++) {
                    data[i].add(new pair(power_values[i][k][m], r_values[i][k][m], k, m,i));
                }
            }
        }

    }

    public static void main(String[] args) throws Exception {
        /*
        String before_prep = "";
        String prep1 = "";
        String prep2 = "";
        String prep3 = "";
        for(int i =4; i < 6; i++) {
            try {
                long sum = 0;
                for (int j = 0; j < 100; j++)
                {
                    Reader rr = new Reader("testfiles-2/test"+String.valueOf(i)+ ".txt");
                rr.remove_impossible_terms();
                rr.remove_IP_dominated();
                rr.remove_LP_dominated();
                    long start = System.nanoTime();
                rr.greedy_LP();
                long time = System.nanoTime() - start;
                sum += time;
                }
                before_prep += String.valueOf((sum/100.0) * Math.pow(10,-6)) + " & ";

            }
            catch (Exception e){
                before_prep += String.valueOf("N.A") + " & ";

            }

        }
        System.out.println(before_prep);
        */


        Reader rr = new Reader("testfiles-2/test4.txt");
        rr.remove_impossible_terms();
        rr.remove_IP_dominated();
        rr.remove_LP_dominated();


        Solution sol = rr.greedy_LP();
        //for(ArrayList l:sol.data)  System.out.println(l);
        double max_r = rr.LP_solver();
        System.out.print("Maximum rate found by LP_solver is : ");System.out.println(max_r);
        System.out.print("Maximum rate found by greedy is : "); System.out.println(sol.Rate);


        System.out.print("Maximum rate found by DP1 is : "); System.out.println(rr.DP_1());

        System.out.print("Maximum rate found by DP2 is : "); System.out.println(rr.DP_2(rr.upper_bnd_Rate));
        System.out.print("Maximum rate found by BB is : ");System.out.println(rr.Braunch_and_bound());


    }
    int total_num_instances(boolean after_lp_filter){
        int c = 0;
        ArrayList<pair>[] data_to_use = data;
        if (after_lp_filter)  data_to_use = LP_filtered_data;
        for (ArrayList dd: data_to_use){
            c += dd.size();
        }
        return c;
    }
    void remove_impossible_terms() throws Exception { //First step of pre-processing
        pair[] mins = new pair[N]; //Array to hold the pairs with minimum power for each channel
        int power_sum_min = 0;
        for (int i = 0; i < N; i++) {
            pair min = Collections.min(data[i], new power_comparator());
            mins[i] = min;
            power_sum_min += min.p;
        }
        if (power_sum_min > this.p) throw new Exception("Impossible problem instance");
        for (int n = 0; n < N; n++) {
            ArrayList<pair> channel_n = data[n];
            HashSet<pair> channel_n_hash = new HashSet<pair>(channel_n);
            for (pair pai : channel_n) {
                if (power_sum_min + pai.p - mins[n].p > this.p) channel_n_hash.remove(p);
            }
            data[n] = new ArrayList(channel_n_hash);
        }

    }

    void remove_IP_dominated() {
        for (int n = 0; n < N; n++) {
            ArrayList<pair> channel_n = data[n];
            HashSet<pair> channel_n_hash = new HashSet<pair>(channel_n); //to remove dominated terms in constant time
            Collections.sort(channel_n, new power_comparator());
            int max_r = channel_n.get(0).r;
            for (int i = 1; i < channel_n.size(); i++) {
                pair curr_pair = channel_n.get(i);
                if (curr_pair.r <= max_r) channel_n_hash.remove(curr_pair);
                max_r = Math.max(max_r, curr_pair.r);
            }
            data[n] = new ArrayList(channel_n_hash);
        }

    }

    boolean is_point_right(pair p1, pair p2, pair p3) { //true if p1 is on the right of the line p3 -> p2

        return (p2.r - p3.r) * (p1.p - p2.p) >= (p1.r - p2.r) * (p2.p - p3.p);


    }

    void remove_LP_dominated() {
        for (int n = 0; n < N; n++) {
            LP_filtered_data[n] = new ArrayList<pair>();
            ArrayList<pair> channel_n = data[n];
            Collections.sort(channel_n, new power_comparator());
            Stack<pair> upper_convex_hull = new Stack<pair>();
            upper_convex_hull.push(channel_n.get(0));
            for (int i = 1; i < channel_n.size(); i++) {
                pair p1 = channel_n.get(i);
                if (upper_convex_hull.size() > 1) {
                    pair p2 = upper_convex_hull.pop();
                    pair p3 = upper_convex_hull.peek();
                    while (!is_point_right(p1, p2, p3) && upper_convex_hull.size() > 1) {
                        p2 = upper_convex_hull.pop();
                        p3 = upper_convex_hull.peek();


                    }
                    if (is_point_right(p1, p2, p3)) upper_convex_hull.push(p2);
                }
                upper_convex_hull.push(p1);

            }
            LP_filtered_data[n] = new ArrayList<pair>(upper_convex_hull);
        }

    }

    void visualize_data(int channel, String additionnal_title, boolean lp_filtred) {
        data_visualiser chart = null;
        if (lp_filtred) {
            chart = new data_visualiser("Instance scatter plot - file :" + filename + additionnal_title, LP_filtered_data, channel);

        }
        else {chart = new data_visualiser("Instance scatter plot - file :" + filename + additionnal_title, data, channel);}
        chart.setSize(800, 400);
        chart.setLocationRelativeTo(null);
        chart.setVisible(true);
    }

    ArrayList<pair> Sort_by_incremental_efficiency() { //returns an array of all pairs apart form the first pair of each channel, sorted by incremental efficiency
        ArrayList<pair> pairs_sorted_eff = new ArrayList<pair>();
        for (int n = 0; n < N; n++) {
            ArrayList<pair> channel_n = LP_filtered_data[n]; //OPERATES ON data after removing LP dominated terms
            Collections.sort(channel_n, new power_comparator());
            pair curr_pair = channel_n.get(0);
            curr_pair.inc_eff = Double.MAX_VALUE;
            curr_pair.inc_power = curr_pair.p;
            curr_pair.inc_rate = curr_pair.r;
            pairs_sorted_eff.add(curr_pair);
            for (int i = 1; i < channel_n.size(); i++) {
                curr_pair = channel_n.get(i);
                pair prev_pair = channel_n.get(i - 1);
                curr_pair.inc_eff = Double.valueOf(curr_pair.r - prev_pair.r) / (curr_pair.p - prev_pair.p);
                curr_pair.inc_rate = (curr_pair.r - prev_pair.r);
                curr_pair.inc_power = (curr_pair.p - prev_pair.p);
                pairs_sorted_eff.add(curr_pair);
            }
        }
        Collections.sort(pairs_sorted_eff,new inc_eff_comparator().reversed());
        return pairs_sorted_eff;
    }

    Solution greedy_LP() {
        ArrayList<pair> sorted_inc = Sort_by_incremental_efficiency();

        int Power_Bud = this.p;
        double Rate = 0;
        ArrayList<sol_pair>[] solutions = new ArrayList[N];
        for (int j = 0; j < N; j++) solutions[j] = new ArrayList<>();
        int i = 0;
        pair curr_pair = null;
        while (i < sorted_inc.size() && sorted_inc.get(i).inc_power <= Power_Bud ) {
            curr_pair = sorted_inc.get(i);
            Power_Bud -= curr_pair.inc_power;
            Rate += curr_pair.inc_rate;
            sol_pair sol = new sol_pair(curr_pair.p,curr_pair.r,curr_pair.user,curr_pair.m,curr_pair.n,1);
            solutions[curr_pair.n].clear();
            solutions[curr_pair.n].add(sol);
            i++;
        }

        if (Power_Bud > 0 && i < sorted_inc.size()) {
            curr_pair = sorted_inc.get(i);
            double x = Double.valueOf(Power_Bud)/curr_pair.inc_power;

            sol_pair sol1 = new sol_pair(curr_pair.p,curr_pair.r,curr_pair.user,curr_pair.m,curr_pair.n,x);
            pair prev_pair = solutions[sol1.n].get(0);
            sol_pair sol2 = new sol_pair(prev_pair.p,prev_pair.r,prev_pair.user,prev_pair.m,prev_pair.n,1-x);
            Rate += x*curr_pair.inc_rate;
            Power_Bud -= x*curr_pair.inc_power;
            solutions[sol1.n].clear();
            solutions[sol1.n].add(sol2);
            solutions[sol1.n].add(sol1);

        }
        return new Solution(Rate,solutions);


    }

    double LP_solver() {
        ArrayList<Double> rates = new ArrayList<Double>();
        ArrayList<Double> powers = new ArrayList<Double>();
        int[] channels_sizes = new int[N];
        for(int i = 0; i < N; i++) {
            for (pair p:data[i]){
                rates.add(Double.valueOf(p.r));
                powers.add(Double.valueOf(p.p));
            }
            if (i==0) channels_sizes[i] = 0;
            else channels_sizes[i] = channels_sizes[i-1] + data[i-1].size();
        }
        LinearProgram lp = new LinearProgram(rates.stream().mapToDouble(Double::doubleValue).toArray());
        lp.addConstraint(new LinearSmallerThanEqualsConstraint(powers.stream().mapToDouble(Double::doubleValue).toArray(), this.p, "power budget constraint"));
        for(int i = 0; i < N; i++) {
            ArrayList<Double> sparse_vec_ch_i = new ArrayList<Double>(); //to hold the list used to build the ith constraint
            Double[] zeros_arr_before = new Double[channels_sizes[i]];
            Arrays.fill(zeros_arr_before,0.0);
            Collections.addAll(sparse_vec_ch_i,zeros_arr_before);
            Double[] ones_arr = new Double[data[i].size()];
            Arrays.fill(ones_arr,1.0);
            Collections.addAll(sparse_vec_ch_i,ones_arr);
            if (i < N-1) {
                Double[] zeros_arr_after = new Double[channels_sizes[N-1] + data[N-1].size() - channels_sizes[i+1]];
                Arrays.fill(zeros_arr_after,0.0);
                Collections.addAll(sparse_vec_ch_i,zeros_arr_after);
            }
            lp.addConstraint(new LinearEqualsConstraint(sparse_vec_ch_i.stream().mapToDouble(Double::doubleValue).toArray(),1.0,"One user per channel constraint" + String.valueOf(i)));

        }
        double[] is_int = new double[rates.size()];
        for (int i = 0; i < rates.size(); i++) {
            is_int[i] = 0;
        }
        lp.setLowerbound(is_int);

        LinearProgramSolver solver  = SolverFactory.newDefault();
        double[] sol = solver.solve(lp);


        double maxrate = 0;
        for (int i = 0; i < rates.size();i++){
            maxrate += sol[i]*rates.get(i);
        }
        return maxrate;

    }

    int DP_1(){ //First Implementation of Dynamic programming
        /*
        L[p-1] holds maximum rate with power budget p
          */
        int[] L = new int[p];
        for(int i = 1; i <= p; i++) {
            int maxr = 0;
            for(pair curr_pair:data[0]){
                if (curr_pair.p <= i) maxr = Math.max(maxr,curr_pair.r);
            }
            L[i-1] = maxr;
        }
        for(int i =1; i < N; i++){
            int[] tempL = new int[p];
            for (int power = 1; power <= p; power++) {
                ArrayList<pair> curr_channel = data[i];
                int max_r = 0;
                for(pair pp:curr_channel) {
                    if (pp.p < power && L[power-pp.p - 1] > 0) max_r = Math.max(max_r,L[power-pp.p - 1] + pp.r);
                }
                tempL[power-1] = max_r;
            }
            L = tempL;
        }
        return L[p-1];
    }

    int DP_2(int U){ //First Implementation of Dynamic programming given upper bound U of rates
        int[] L = new int[U];
        for(int i = 1; i <=U; i++) {
            int minp = 0;
            for(pair curr_pair:data[0]){
                if (curr_pair.r == i){
                    if(minp ==0) minp = curr_pair.p;
                    else minp = Math.min(minp,curr_pair.p);
                }
            }
            L[i-1] = minp;
        }
        for(int n =1; n < N; n++){
            int[] tempL = new int[U];
            for (int rate = 1; rate <=U; rate++) {
                ArrayList<pair> curr_channel = data[n];
                int min_p = 0;
                for(pair pp:curr_channel) {
                    if (pp.r < rate && L[rate-pp.r - 1] > 0){
                        if (min_p == 0) min_p = L[rate-pp.r - 1] + pp.p;
                        else min_p = Math.min(min_p,L[rate-pp.r - 1] + pp.p);
                    }
                }
                tempL[rate-1] = min_p;
            }
            L = tempL;
        }
        for(int r = U; r >= 1; r--) {
            if (L[r-1] <= p && L[r-1] > 0) return r;
        }
        return -1;
    }

    Bounds Greedy_bound(int curr_channel,int Power_Bud , ArrayList<pair> sorted_inc) {
        double Rate = 0;
        int i = 0;
        pair curr_pair = null;
        while (i < sorted_inc.size() && sorted_inc.get(i).inc_power <= Power_Bud ) {
            curr_pair = sorted_inc.get(i);
            if (curr_pair.n >= curr_channel) {
                Power_Bud -= curr_pair.inc_power;
                Rate += curr_pair.inc_rate;
                sol_pair sol = new sol_pair(curr_pair.p, curr_pair.r, curr_pair.user, curr_pair.m, curr_pair.n, 1);
            }
            i++;
        }
        int LB = (int) Rate;
        if (Power_Bud > 0 && i < sorted_inc.size()) {
            curr_pair = sorted_inc.get(i);
            double x = Double.valueOf(Power_Bud)/curr_pair.inc_power;
            Rate += x*curr_pair.inc_rate;
            Power_Bud -= x*curr_pair.inc_power;

        }
        return new Bounds(Rate,LB);
    }
    int Braunch_and_bound(){
        Bounds curbound = new Bounds(Double.MAX_VALUE,Integer.MIN_VALUE);
        BB(0,0,0,curbound, Sort_by_incremental_efficiency());
        return curbound.LB;
    }
    void BB(int curr_channel, int Power_used, int rate_achieved, Bounds curr_bounds, ArrayList<pair> sorted_inc) {
        if (Power_used >= this.p) return;
        Bounds braunch_bound = Greedy_bound(curr_channel,this.p - Power_used,sorted_inc);
        /*
        System.out.print(braunch_bound.UB + rate_achieved);
        System.out.print(",");
        System.out.print(braunch_bound.LB + rate_achieved);
        System.out.print(",");
        System.out.println(curr_channel);
*/
        if (braunch_bound.UB + rate_achieved <= curr_bounds.LB) return;
        curr_bounds.LB = Math.max(curr_bounds.LB,braunch_bound.LB + rate_achieved);
        curr_bounds.UB = Math.min(curr_bounds.UB,braunch_bound.UB + rate_achieved);
        for (pair pai:data[curr_channel]){
            if (pai.p + Power_used > this.p) continue;
            if (curr_channel < N-1) BB(curr_channel + 1, Power_used + pai.p, rate_achieved + pai.r,curr_bounds,sorted_inc);
            else curr_bounds.LB = Math.max(curr_bounds.LB,rate_achieved + pai.r);
        }
    }


}
class power_comparator implements Comparator{
    public int compare(Object o1, Object o2) {
        pair pair1 = (pair) o1;
        pair pair2 = (pair) o2;
        if (pair1.p > pair2.p) return 1;
        else if(pair1.p < pair2.p) return -1;
        else {
            if (pair1.r > pair2.r) return -1;
            else if (pair1.r < pair2.r) return 1;
            else return 0;
        }

    }
}
class inc_eff_comparator implements  Comparator {
    public int compare(Object o1, Object o2) {
        pair pair1 = (pair) o1;
        pair pair2 = (pair) o2;
        return Double.compare(pair1.inc_eff,pair2.inc_eff);
    }
}
class Bounds {
    double UB;
    int LB;
    public Bounds(double ub, int lb) {
        UB = ub;
        LB = lb;
    }
}



