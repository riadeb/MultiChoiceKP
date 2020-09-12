import containers.Solution;

import java.util.Scanner;

public class Main {
    public static void main(String[] args) throws Exception {
        Scanner in = new Scanner(System.in);
        System.out.println("Pleaser enter the file path of your instance, or X to exit");
        String s = in.nextLine();

        while (!s.equals("X")){


            String file = s;
            MCKP_Solver rr = new MCKP_Solver(file);
            try {
                rr.remove_impossible_terms();
            }
            catch (Exception e){
                System.out.println("Impossible instance");
                System.out.println("Pleaser enter the file path of your instance, or X to exit");
                s = in.nextLine();
                continue;
            }

            rr.remove_IP_dominated();
            rr.visualize_data(0," After IP_dom removed",false);
            rr.remove_LP_dominated();
            rr.visualize_data(0," After LP_dom removed",true);
            Solution sol = rr.greedy_LP();
            //for(ArrayList l:sol.data)  System.out.println(l);
            double max_r = rr.LP_solver();
            System.out.print("Maximum rate found by LP_solver is : ");System.out.println(max_r);
            System.out.print("Maximum rate found by greedy is : "); System.out.println(sol.Rate);
            System.out.print("Maximum rate found by DP1 is : "); System.out.println(rr.DP_1());
            System.out.print("Maximum rate found by DP2 is : "); System.out.println(rr.DP_2(rr.upper_bnd_Rate));
            System.out.print("Maximum rate found by BB is : ");System.out.println(rr.Braunch_and_bound());

            System.out.println("Pleaser enter the file path of your instance, or X to exit");
            s = in.nextLine();
        }


    }

}
