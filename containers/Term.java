package containers;

import java.util.ArrayList;

public class Term {
    public int p;
    public int r;
    public int n; //channel index
    public double inc_eff; //equals inc_rate/inc_eff
    public int inc_power;
    public int inc_rate;
    public double weight;
    public  Term(int p, int r, int i) {
        this.p = p;
        this.r= r;
        this.n = i;
        this.weight = 1;
    }
    public  Term(int p, int r, int i, double w) {
        this.p = p;
        this.r= r;
        this.n = i;
        this.weight = w;
    }


    //abstract public boolean equals(Object obj);
    //abstract public String toString();
}