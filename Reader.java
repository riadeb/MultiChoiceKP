import java.util.*;
import  java.io.File;


public class Reader {
    int p;
    int M;
    int K;
    int N;
    HashSet [] data; //array of hashsets (to remove dominated elements in constant time)
    public Reader(String filename) throws Exception {
        File file = new File(filename);
        Scanner sc = new Scanner(file);
        this.N = Double.valueOf(sc.nextLine()).intValue();
        this.M = Double.valueOf(sc.nextLine()).intValue();
        this.K = Double.valueOf(sc.nextLine()).intValue();
        this.p = Double.valueOf(sc.nextLine()).intValue();
        data = new HashSet[N];

        int[][][] power_values = new int[N][K][M];
        int[][][] r_values = new int[N][K][M];

        for (int i = 0; i < N; i++) {
            data[i] = new HashSet<pair>();
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
        Reader rr = new Reader("testfiles-2/test5.txt");
        System.out.println(rr.p);
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
}
