import java.util.*;
import  java.io.File;



public class Reader {
    int p;
    int M;
    int K;
    int N;
    String filename;
    ArrayList<pair> [] data;
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
            for(int k =0; k < K; k++) {
                String[] newline = sc.nextLine().stripLeading().split("   ");
                for(int m = 0; m < M; m++){

                        power_values[i][k][m] = Double.valueOf(newline[m]).intValue();

                }
            }
        }
        for (int i = 0; i < N; i++) {
            for(int k =0; k < K; k++) {
                String[] newline = sc.nextLine().stripLeading().split("   ");
                for(int m = 0; m < M; m++){
                    r_values[i][k][m] = Double.valueOf(newline[m]).intValue();
                }
            }
        }
        for (int i = 0; i < N; i++) {
            for(int k =0; k < K; k++) {
                for(int m = 0; m < M; m++){
                    data[i].add(new pair(power_values[i][k][m],r_values[i][k][m],k,m));
                }
            }
        }

    }
    public static void main(String[] args) throws Exception {
        Reader rr = new Reader("testfiles-2/test4.txt");
        rr.visualize_data(28," Before preprocessing");
        rr.remove_impossible_terms();
        rr.visualize_data(28," After removing impossible terms");
        rr.remove_IP_dominated();
        rr.visualize_data(28," After removing IP dominated terms");
        rr.remove_LP_dominated();
        rr.visualize_data(28," After removing LP dominated terms");

    }
    void remove_impossible_terms() throws Exception{ //First step of pre-processing
        pair[] mins = new pair[N]; //Array to hold the pairs with minimum power for each channel
        int power_sum_min = 0;
        for(int i = 0; i < N; i++){
            pair min = Collections.min(data[i],new power_comparator());
            mins[i] = min;
            power_sum_min += min.p;
        }
        if (power_sum_min > this.p) throw new Exception("Impossible problem instance");
        for (int n = 0; n < N; n++) {
            ArrayList<pair> channel_n = data[n];
            HashSet<pair> channel_n_hash = new HashSet<pair>(channel_n);
            for (pair pai:channel_n) {
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
    boolean is_point_right(pair p1,pair p2,pair p3){ //true if p1 is on the right of the line p3 -> p2

            return (p2.r - p3.r)*(p1.p - p2.p)>= (p1.r - p2.r)*(p2.p - p3.p) ;


    }
    void remove_LP_dominated(){
        for (int n = 0; n < N; n++) {
            ArrayList<pair> channel_n = data[n];
            Collections.sort(channel_n ,new power_comparator());
            Stack<pair> upper_convex_hull = new Stack<pair>();
            upper_convex_hull.push(channel_n.get(0));
            for(int i=1; i < channel_n.size(); i++){
                pair p1 = channel_n.get(i);
                if (upper_convex_hull.size() > 1) {
                    pair p2 = upper_convex_hull.pop();
                    pair p3 = upper_convex_hull.peek();
                    while(!is_point_right(p1,p2,p3) && upper_convex_hull.size() > 1) {
                        p2 = upper_convex_hull.pop();
                        p3 = upper_convex_hull.peek();


                    }
                    if (is_point_right(p1,p2,p3)) upper_convex_hull.push(p2);
                }
                upper_convex_hull.push(p1);

            }
            data[n] = new ArrayList<pair>(upper_convex_hull);
        }
    }
    void visualize_data(int channel,String additionnal_title){
       data_visualiser chart = new data_visualiser("Instance scatter plot - file :" + filename + additionnal_title,data,channel);
        chart.setSize(800, 400);
        chart.setLocationRelativeTo(null);
        chart.setVisible(true);
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

class pair {
    int p;
    int r;
    int user; // user index
    int m; //power index

    public  pair(int p, int r, int u, int m) {
        this.p = p;
        this.r= r;
        this.user = u;
        this.m = m;
    }

    @Override
    public boolean equals(Object obj) {
        pair newp = (pair) obj;
        return (user == newp.user && m == newp.m);

    }
    @Override
    public String toString(){
        return "("+String.valueOf(this.p) + ","+String.valueOf(this.r)+")";
    }
}
