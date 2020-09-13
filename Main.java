import containers.Solution;
import containers.Term;
import solver.MCKP_Solver;

import java.io.File;
import java.util.LinkedList;
import java.util.Scanner;

public class Main {
    public static void main(String[] args) throws Exception {
        if(args.length < 1) {
            System.out.println("Please enter a valid argument");
            System.exit(0);
        }
        if(args.length == 2 && args[0].equals("-r")){
            String file = args[1];
            runtime_comparator(file);
        }
        else if (args.length == 1) {
            String file = args[0];
            MCKP_Solver rr = read_data(file);

                rr.remove_impossible_terms();

            rr.remove_IP_dominated();
            rr.visualize_data(0, " After IP_dom removed", false);
            rr.remove_LP_dominated();
            rr.visualize_data(0, " After LP_dom removed", true);
            Solution sol = rr.greedy_LP();
            //for(ArrayList l:sol.data)  System.out.println(l);
            double max_r = rr.LP_solver();
            System.out.print("Maximum rate found by LP_solver is : ");
            System.out.println(max_r);
            System.out.print("Maximum rate found by greedy is : ");
            System.out.println(sol.profit);
            System.out.print("Maximum rate found by DP1 is : ");
            System.out.println(rr.DP_1());
            System.out.print("Maximum rate found by DP2 is : ");
            System.out.println(rr.DP_2(rr.upper_bnd_profit()));
            System.out.print("Maximum rate found by BB is : ");
            System.out.println(rr.Braunch_and_bound());

        }
        else System.out.println("Invalid argument");
    }
     static MCKP_Solver read_data(String filename) throws Exception {
        File file = new File(filename);
        Scanner sc = new Scanner(file);
        int N = Double.valueOf(sc.nextLine()).intValue();
        int M = Double.valueOf(sc.nextLine()).intValue();
        int K = Double.valueOf(sc.nextLine()).intValue();
        int p = Double.valueOf(sc.nextLine()).intValue();
        LinkedList<Term> [] data = new LinkedList[N];


        int[][][] power_values = new int[N][K][M];
        int[][][] r_values = new int[N][K][M];

        for (int i = 0; i < N; i++) {
            data[i] = new LinkedList<Term>();
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
        }
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < K; k++) {
                for (int m = 0; m < M; m++) {
                    data[i].add(new Term(power_values[i][k][m], r_values[i][k][m],i));
                }
            }
        }

        return new MCKP_Solver(p, N,data);

    }
    static void runtime_comparator(String filename) throws Exception { //used to make testes and calculate runtime

        String[] algos = new String[]{"dp1","dp2","bb"};
        System.out.println("Algorithm | Runtime (ms) | Result");
        for(String alg:algos) {
            String times = "";
            String rates = "";
            String power = "";
            long sum = 0;
            int DP1 = 0;
            int power_used = 0;
            MCKP_Solver rr = read_data(filename);
            try {
                rr.remove_impossible_terms();
            } catch (Exception e) {
                times += String.valueOf("N.A") + " & ";
                rates += String.valueOf("N.A") + " & ";
                continue;

            }
            rr.remove_IP_dominated();
            rr.remove_LP_dominated();
            for (int j = 0; j < 100; j++) {
                long start = System.nanoTime();
                switch (alg){
                    case "dp1":
                        DP1 = rr.DP_1();
                        break;
                    case "dp2":
                        DP1 = rr.DP_2(rr.upper_bnd_profit());
                        break;
                    case "bb":
                         DP1 = rr.Braunch_and_bound();

                }
                long time = System.nanoTime() - start;
                sum += time;
            }
            times += String.valueOf((sum / 100.0) * Math.pow(10, -6)) + " | ";
            rates += String.valueOf(DP1) ;


        System.out.print(alg + " | ");
        System.out.print(times);
        System.out.println(rates);
        }

    }



}
