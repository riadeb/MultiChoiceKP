package solver;

import java.util.*;

import plotter.Plot;
import containers.Solution;
import containers.Term;
import scpsolver.lpsolver.LinearProgramSolver;
import scpsolver.lpsolver.SolverFactory;
import scpsolver.problems.LinearProgram;
import  scpsolver.constraints.*;


public class MCKP_Solver {
    int c;
    int m;
    LinkedList<Term>[] data;
    LinkedList<Term>[] LP_filtered_data;

    public MCKP_Solver(int _c, int _m, LinkedList<Term>[] _data) throws Exception {
        data = _data;
        this.c = _c;
        this.m = _m;
        LP_filtered_data = new LinkedList[m];


    }

    int total_num_instances(boolean after_lp_filter) { //returns total number of instances
        int res = 0;
        LinkedList<Term>[] data_to_use = data;
        if (after_lp_filter) data_to_use = LP_filtered_data;
        for (LinkedList dd : data_to_use) {
            res += dd.size();
        }
        return res;
    }

    public int upper_bnd_profit() { // returns an upper bound for the maximum profit achievable, by taking the biggest one in every _class
        int res = 0;
        for (int i = 0; i < m; i++) {
            int maxr = 0;
            for (Term t : data[i]) {
                maxr = Math.max(maxr, t.profit);
            }

            res += maxr;
        }
        return res;
    }


    public void remove_impossible_terms() throws Exception { //First step of pre-processing
        Term[] mins = new Term[m]; //Array to hold the pairs with minimum weight for each _class
        int weight_sum_min = 0;
        for (int i = 0; i < m; i++) {
            Term min = Collections.min(data[i], new weight_comparator());
            mins[i] = min;
            weight_sum_min += min.weight;
        }
        if (weight_sum_min > this.c) {


            throw new Exception("Impossible problem instance"); //if minimum weight in each _class exceeds budget, then problem is impossible
        }
        for (int n = 0; n < m; n++) {

            LinkedList<Term> _class_n = data[n];
            final int offset = weight_sum_min - mins[n].weight;
            _class_n.removeIf(t -> offset + t.weight > this.c);
        }

    }

    public void remove_IP_dominated() {
        for (int n = 0; n < m; n++) {
            LinkedList<Term> _class_n = data[n];
            Collections.sort(_class_n, new weight_comparator());
            int max_r = _class_n.get(0).profit;
            Iterator<Term> it = _class_n.iterator();
            if(!it.hasNext()) continue;
            else it.next();
            for (; it.hasNext(); ) {
                Term curr = it.next();
                if (curr.profit <= max_r) it.remove();
                max_r = Math.max(max_r, curr.profit);
            }
        }

    }

    boolean is_point_right(Term p1, Term p2, Term p3) { //true if p1 is on the right of the line p3 -> p2

        return (p2.profit - p3.profit) * (p1.weight - p2.weight) >= (p1.profit - p2.profit) * (p2.weight - p3.weight);


    }

    public void remove_LP_dominated() {
        for (int n = 0; n < m; n++) {
            //  LP_filtered_data[n] = new ArrayList<Term>();
            LinkedList<Term> _class_n = data[n];
            Collections.sort(_class_n, new weight_comparator());
            Stack<Term> upper_convex_hull = new Stack<Term>();
            upper_convex_hull.push(_class_n.get(0));
            for (int i = 1; i < _class_n.size(); i++) {
                Term p1 = _class_n.get(i);
                if (upper_convex_hull.size() > 1) {
                    Term p2 = upper_convex_hull.pop();
                    Term p3 = upper_convex_hull.peek();
                    while (!is_point_right(p1, p2, p3) && upper_convex_hull.size() > 1) {
                        p2 = upper_convex_hull.pop();
                        p3 = upper_convex_hull.peek();


                    }
                    if (is_point_right(p1, p2, p3)) upper_convex_hull.push(p2);
                }
                upper_convex_hull.push(p1);

            }
            LP_filtered_data[n] = new LinkedList<Term>(upper_convex_hull);
        }

    }

    public void visualize_data(int _class, String additionnal_title, boolean lp_filtred) {
        Plot chart = null;
        if (lp_filtred) {
            chart = new Plot("Instance scatter plot - " + additionnal_title, LP_filtered_data, _class);

        } else {
            chart = new Plot("Instance scatter plot - "  + additionnal_title, data, _class);
        }
        chart.setSize(800, 400);
        chart.setLocationRelativeTo(null);
        chart.setVisible(true);
    }

    LinkedList<Term> Sort_by_incremental_efficiency() { //returns an array of all pairs apart form the first pair of each _class, sorted by incremental efficiency
        LinkedList<Term> pairs_sorted_eff = new LinkedList<Term>();
        for (int n = 0; n < m; n++) {
            LinkedList<Term> _class_n = LP_filtered_data[n]; //OPEprofitS ON data after removing LP dominated terms
            Collections.sort(_class_n, new weight_comparator());
            Term curr_pair = _class_n.get(0);
            curr_pair.inc_eff = Double.MAX_VALUE;
            curr_pair.inc_weight = curr_pair.weight;
            curr_pair.inc_profit = curr_pair.profit;
            pairs_sorted_eff.add(curr_pair);
            for (int i = 1; i < _class_n.size(); i++) {
                curr_pair = _class_n.get(i);
                Term prev_pair = _class_n.get(i - 1);
                curr_pair.inc_eff = Double.valueOf(curr_pair.profit - prev_pair.profit) / (curr_pair.weight - prev_pair.weight);
                curr_pair.inc_profit = (curr_pair.profit - prev_pair.profit);
                curr_pair.inc_weight = (curr_pair.weight - prev_pair.weight);
                pairs_sorted_eff.add(curr_pair);
            }
        }
        Collections.sort(pairs_sorted_eff, new inc_eff_comparator().reversed());
        return pairs_sorted_eff;
    }

    public Solution greedy_LP() {
        LinkedList<Term> sorted_inc = Sort_by_incremental_efficiency();

        int weight_Bud = this.c;
        double profit = 0;
        ArrayList<Term>[] solutions = new ArrayList[m];
        for (int j = 0; j < m; j++) solutions[j] = new ArrayList<>();
        int i = 0;
        Term curr_pair = null;
        while (i < sorted_inc.size() && sorted_inc.get(i).inc_weight <= weight_Bud) {
            curr_pair = sorted_inc.get(i);
            weight_Bud -= curr_pair.inc_weight;
            profit += curr_pair.inc_profit;
            Term sol = new Term(curr_pair.weight, curr_pair.profit, curr_pair.n, 1);
            solutions[curr_pair.n].clear();
            solutions[curr_pair.n].add(sol);
            i++;
        }

        if (weight_Bud > 0 && i < sorted_inc.size()) {
            curr_pair = sorted_inc.get(i);
            double x = Double.valueOf(weight_Bud) / curr_pair.inc_weight;

            Term sol1 = new Term(curr_pair.weight, curr_pair.profit, curr_pair.n, x);
            Term prev_pair = solutions[sol1.n].get(0);
            Term sol2 = new Term(prev_pair.weight, prev_pair.profit, prev_pair.n, 1 - x);
            profit += x * curr_pair.inc_profit;
            weight_Bud -= x * curr_pair.inc_weight;
            solutions[sol1.n].clear();
            solutions[sol1.n].add(sol2);
            solutions[sol1.n].add(sol1);

        }
        return new Solution(profit, solutions);


    }

    public double LP_solver() {
        ArrayList<Double> profits = new ArrayList<Double>();
        ArrayList<Double> weights = new ArrayList<Double>();
        int[] _classs_sizes = new int[m];
        for (int i = 0; i < m; i++) {
            for (Term p : data[i]) {
                profits.add(Double.valueOf(p.profit));
                weights.add(Double.valueOf(p.weight));
            }
            if (i == 0) _classs_sizes[i] = 0;
            else _classs_sizes[i] = _classs_sizes[i - 1] + data[i - 1].size();
        }
        LinearProgram lp = new LinearProgram(profits.stream().mapToDouble(Double::doubleValue).toArray());
        lp.addConstraint(new LinearSmallerThanEqualsConstraint(weights.stream().mapToDouble(Double::doubleValue).toArray(), this.c, "weight budget constraint"));
        for (int i = 0; i < m; i++) {
            ArrayList<Double> sparse_vec_ch_i = new ArrayList<Double>(); //to hold the list used to build the ith constraint
            Double[] zeros_arr_before = new Double[_classs_sizes[i]];
            Arrays.fill(zeros_arr_before, 0.0);
            Collections.addAll(sparse_vec_ch_i, zeros_arr_before);
            Double[] ones_arr = new Double[data[i].size()];
            Arrays.fill(ones_arr, 1.0);
            Collections.addAll(sparse_vec_ch_i, ones_arr);
            if (i < m - 1) {
                Double[] zeros_arr_after = new Double[_classs_sizes[m - 1] + data[m - 1].size() - _classs_sizes[i + 1]];
                Arrays.fill(zeros_arr_after, 0.0);
                Collections.addAll(sparse_vec_ch_i, zeros_arr_after);
            }
            lp.addConstraint(new LinearEqualsConstraint(sparse_vec_ch_i.stream().mapToDouble(Double::doubleValue).toArray(), 1.0, "One user per _class constraint" + String.valueOf(i)));

        }
        double[] is_int = new double[profits.size()];
        for (int i = 0; i < profits.size(); i++) {
            is_int[i] = 0;
        }
        lp.setLowerbound(is_int);

        LinearProgramSolver solver = SolverFactory.newDefault();
        double[] sol = solver.solve(lp);


        double maxprofit = 0;
        for (int i = 0; i < profits.size(); i++) {
            maxprofit += sol[i] * profits.get(i);
        }
        return maxprofit;

    }

    public int DP_1() { //First Implementation of Dynamic programming
        /*
        L[p-1] holds maximum profit with weight budget p
          */
        int[] L = new int[c];
        for (int i = 1; i <= c; i++) {
            int maxr = 0;
            for (Term curr_pair : data[0]) {
                if (curr_pair.weight <= i) maxr = Math.max(maxr, curr_pair.profit);
            }
            L[i - 1] = maxr;
        }
        for (int i = 1; i < m; i++) {
            int[] tempL = new int[c];
            for (int weight = 1; weight <= c; weight++) {
                LinkedList<Term> curr__class = data[i];
                int max_r = 0;
                for (Term pp : curr__class) {
                    if (pp.weight < weight && L[weight - pp.weight - 1] > 0) max_r = Math.max(max_r, L[weight - pp.weight - 1] + pp.profit);
                }
                tempL[weight - 1] = max_r;
            }
            L = tempL;
        }
        return L[c - 1];
    }

    public int DP_2(int U) { //First Implementation of Dynamic programming given upper bound U of profits
        int[] L = new int[U];
        for (int i = 1; i <= U; i++) {
            int minp = 0;
            for (Term curr_pair : data[0]) {
                if (curr_pair.profit == i) {
                    if (minp == 0) minp = curr_pair.weight;
                    else minp = Math.min(minp, curr_pair.weight);
                }
            }
            L[i - 1] = minp;
        }
        for (int n = 1; n < m; n++) {
            int[] tempL = new int[U];
            for (int profit = 1; profit <= U; profit++) {
                LinkedList<Term> curr__class = data[n];
                int min_p = 0;
                for (Term pp : curr__class) {
                    if (pp.profit < profit && L[profit - pp.profit - 1] > 0) {
                        if (min_p == 0) min_p = L[profit - pp.profit - 1] + pp.weight;
                        else min_p = Math.min(min_p, L[profit - pp.profit - 1] + pp.weight);
                    }
                }
                tempL[profit - 1] = min_p;
            }
            L = tempL;
        }
        for (int r = U; r >= 1; r--) {
            if (L[r - 1] <= c && L[r - 1] > 0) return r;
        }
        return -1;
    }

    Bounds Greedy_bound(int curr__class, int weight_Bud, LinkedList<Term> sorted_inc) {
        double profit = 0;
        int i = 0;
        Term curr_pair = null;
        while (i < sorted_inc.size() && sorted_inc.get(i).inc_weight <= weight_Bud) {
            curr_pair = sorted_inc.get(i);
            if (curr_pair.n >= curr__class) {
                weight_Bud -= curr_pair.inc_weight;
                profit += curr_pair.inc_profit;
                Term sol = new Term(curr_pair.weight, curr_pair.profit, curr_pair.n, 1);
            }
            i++;
        }
        int LB = (int) profit;
        if (weight_Bud > 0 && i < sorted_inc.size()) {
            curr_pair = sorted_inc.get(i);
            double x = Double.valueOf(weight_Bud) / curr_pair.inc_weight;
            profit += x * curr_pair.inc_profit;
            weight_Bud -= x * curr_pair.inc_weight;

        }
        return new Bounds(profit, LB);
    }

    public int Braunch_and_bound() {
        LinkedList<Term> so = Sort_by_incremental_efficiency();
        Bounds curbound = Greedy_bound(0, this.c, so);
        BB(0, 0, 0, curbound, so);
        return curbound.LB;
    }

    public void BB(int curr__class, int weight_used, int profit_achieved, Bounds curr_bounds, LinkedList<Term> sorted_inc) {
        if (weight_used >= this.c) return;
        Bounds braunch_bound;
        /*
        System.out.print(braunch_bound.UB + profit_achieved);
        System.out.print(",");
        System.out.print(braunch_bound.LB + profit_achieved);
        System.out.print(",");
        System.out.println(curr__class);
*/

        for (Term pai : data[curr__class]) {
            if (pai.weight + weight_used > this.c) continue;
            if (curr__class < m - 1) {
                braunch_bound = Greedy_bound(curr__class + 1, this.c - weight_used - pai.weight, sorted_inc);
                if (braunch_bound.UB + profit_achieved + pai.profit > curr_bounds.LB) {
                    curr_bounds.LB = Math.max(curr_bounds.LB, braunch_bound.LB + profit_achieved + pai.profit);
                    BB(curr__class + 1, weight_used + pai.weight, profit_achieved + pai.profit, curr_bounds, sorted_inc);
                }
            } else curr_bounds.LB = Math.max(curr_bounds.LB, profit_achieved + pai.profit);
        }
    }

    public int BB_with_queue() {
        LinkedList<Term> sorted_inc = Sort_by_incremental_efficiency();
        Queue<par> dd = new LinkedList<>();
        dd.add(new par(0, 0, 0));
        Bounds curr_bounds = Greedy_bound(0, 0, sorted_inc);
        while (dd.size() > 0) {

            par c = dd.remove();
            for (Term pai : data[c.curr__class]) {
                if (pai.weight + c.weight_used > this.c) continue;
                if (c.curr__class < m - 1) {
                    Bounds braunch_bound = Greedy_bound(c.curr__class + 1, this.c - c.weight_used - pai.weight, sorted_inc);
                    if (braunch_bound.UB + c.profit_achieved + pai.profit > curr_bounds.LB) {
                        dd.add(new par(c.curr__class + 1, c.weight_used + pai.weight, c.profit_achieved + pai.profit));
                    }
                    curr_bounds.LB = Math.max(curr_bounds.LB, c.profit_achieved + pai.profit + braunch_bound.LB);
                } else curr_bounds.LB = Math.max(curr_bounds.LB, c.profit_achieved + pai.profit);
            }

        }
        return curr_bounds.LB;
    }
}

class par {
    int curr__class;
    int weight_used;
    int profit_achieved;
    public par(int c, int p, int r) {
        curr__class = c;
        weight_used = p;
        profit_achieved = r;
    }
} //class used to hold nodes parameters
class weight_comparator implements Comparator{
    public int compare(Object o1, Object o2) {
        Term pair1 = (Term) o1;
        Term pair2 = (Term) o2;
        if (pair1.weight > pair2.weight) return 1;
        else if(pair1.weight < pair2.weight) return -1;
        else {
            if (pair1.profit > pair2.profit) return -1;
            else if (pair1.profit < pair2.profit) return 1;
            else return 0;
        }

    }
}
class inc_eff_comparator implements  Comparator {
    public int compare(Object o1, Object o2) {
        Term pair1 = (Term) o1;
        Term pair2 = (Term) o2;
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
} //used to hold upper and lower bounds for problem



