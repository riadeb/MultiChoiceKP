import java.util.*;
import  java.io.File;


public class Reader {
    int p;
    int M;
    int K;
    int N;
    ArrayList [] data;
    public Reader(String filename) throws Exception {
        File file = new File(filename);
        Scanner sc = new Scanner(file);
        this.N = Double.valueOf(sc.nextLine()).intValue();
        this.M = Double.valueOf(sc.nextLine()).intValue();
        this.K = Double.valueOf(sc.nextLine()).intValue();
        this.p = Double.valueOf(sc.nextLine()).intValue();
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
        rr.remmove_IP_dominated();
        System.out.println(rr.data[0]);
    }
    void remmove_IP_dominated(){
        for (int n = 0; n < N; n++) {
            ArrayList<pair> channel_n = data[n];
            HashSet<pair> channel_n_hash = new HashSet<pair>(channel_n); //to remove dominated terms in constant time
            Collections.sort(channel_n ,new power_comparator());
            int max_r = channel_n.get(0).r;
            for(int i = 1; i < channel_n.size(); i++) {
                pair curr_pair = channel_n.get(i);
                if (curr_pair.r <= max_r) channel_n_hash.remove(curr_pair);
                max_r = Math.max(max_r,curr_pair.r);
            }
            data[n] = new ArrayList(channel_n_hash);
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
            if (pair1.r > pair2.r) return 1;
            else if (pair1.r < pair2.r) return -1;
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
}
