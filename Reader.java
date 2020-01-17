import java.lang.reflect.Array;
import java.util.*;
import  java.io.File;



public class Reader {
    int p;
    int M;
    int K;
    int N;
    String filename;
    ArrayList<pair>[] data;

    public Reader(String filename) throws Exception {
        File file = new File(filename);
        Scanner sc = new Scanner(file);
        this.N = Double.valueOf(sc.nextLine()).intValue();
        this.M = Double.valueOf(sc.nextLine()).intValue();
        this.K = Double.valueOf(sc.nextLine()).intValue();
        this.p = Double.valueOf(sc.nextLine()).intValue();
        this.filename = filename;
        data = new ArrayList[N];

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
            for (int k = 0; k < K; k++) {
                String[] newline = sc.nextLine().stripLeading().split("   ");
                for (int m = 0; m < M; m++) {
                    r_values[i][k][m] = Double.valueOf(newline[m]).intValue();
                }
            }
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
        Reader rr = new Reader("testfiles-2/test1.txt");
        rr.remove_impossible_terms();
        rr.remove_IP_dominated();
        rr.remove_LP_dominated();
        rr.visualize_data(3, " After removing LP dominated terms");
        Solution sol = rr.greedy_LP();
        System.out.print("Maximum rate is : "); System.out.println(sol.Rate);
        for(ArrayList l:sol.data)  System.out.println(l);

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
            data[n] = new ArrayList<pair>(upper_convex_hull);
        }
    }

    void visualize_data(int channel, String additionnal_title) {
        data_visualiser chart = new data_visualiser("Instance scatter plot - file :" + filename + additionnal_title, data, channel);
        chart.setSize(800, 400);
        chart.setLocationRelativeTo(null);
        chart.setVisible(true);
    }

    ArrayList<pair> Sort_by_incremental_efficiency() { //returns an array of all pairs apart form the first pair of each channel, sorted by incremental efficiency
        ArrayList<pair> pairs_sorted_eff = new ArrayList<pair>();
        for (int n = 0; n < N; n++) {
            ArrayList<pair> channel_n = data[n];
            Collections.sort(channel_n, new power_comparator());
            pair curr_pair = channel_n.get(0);
            curr_pair.inc_eff = Double.MAX_VALUE;
            curr_pair.inc_power = curr_pair.p;
            curr_pair.inc_rate = curr_pair.r;
            pairs_sorted_eff.add(curr_pair);
            for (int i = 1; i < channel_n.size(); i++) {
                curr_pair = channel_n.get(i);
                pair prev_pair = channel_n.get(i - 1);
                curr_pair.inc_eff = (curr_pair.r - prev_pair.r) / (curr_pair.p - prev_pair.p);
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
        while (i < sorted_inc.size() && sorted_inc.get(i).p <= Power_Bud ) {
            curr_pair = sorted_inc.get(i);
            Power_Bud -= curr_pair.inc_power;
            Rate += curr_pair.inc_rate;
            sol_pair sol = new sol_pair(curr_pair.p,curr_pair.r,curr_pair.user,curr_pair.m,curr_pair.n,1);
            solutions[curr_pair.n].clear();
            solutions[curr_pair.n].add(sol);
            i++;
        }
        if (Power_Bud > 0 && i < sorted_inc.size()) {
            double x = Power_Bud/curr_pair.inc_power;

            sol_pair sol1 = new sol_pair(curr_pair.p,curr_pair.r,curr_pair.user,curr_pair.m,curr_pair.n,x);
            pair prev_pair = solutions[sol1.n].get(0);
            sol_pair sol2 = new sol_pair(prev_pair.p,prev_pair.r,prev_pair.user,prev_pair.m,prev_pair.n,1-x);
            Rate += x*curr_pair.inc_rate;
            Power_Bud -= x*curr_pair.inc_power;
            solutions[sol1.n].clear();
            solutions[sol1.n].add(sol2);
            solutions[sol1.n].add(sol1);
            System.out.println(Power_Bud);

        }
        return new Solution(Rate,solutions);


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



