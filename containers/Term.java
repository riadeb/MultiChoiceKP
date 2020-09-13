package containers;

public class Term {
    public int weight;
    public int profit;
    public int n; //class index
    public double inc_eff;
    public int inc_weight;
    public int inc_profit;
    public double fract;
    public  Term(int weight, int profit, int i) {
        this.weight = weight;
        this.profit = profit;
        this.n = i;
        this.fract = 1;
    }
    public  Term(int weight, int profit, int i, double f) {
        this.weight = weight;
        this.profit = profit;
        this.n = i;
        this.fract = f;
    }


    //abstract public boolean equals(Object obj);
    //abstract public String toString();
}