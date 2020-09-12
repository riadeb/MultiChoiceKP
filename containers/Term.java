package containers;

import java.util.ArrayList;

public class Term {
    public int p;
    public int r;
    public int user; // user index
    public int m; //power index
    public int n; //channel index
    public double inc_eff; //equals inc_rate/inc_eff
    public int inc_power;
    public int inc_rate;
    public double weight;
    public  Term(int p, int r, int u, int m, int i) {
        this.p = p;
        this.r= r;
        this.user = u;
        this.m = m;
        this.n = i;
        this.weight = 1;
    }
    public  Term(int p, int r, int u, int m, int i, double w) {
        this.p = p;
        this.r= r;
        this.user = u;
        this.m = m;
        this.n = i;
        this.weight = w;
    }


    @Override
    public boolean equals(Object obj) {
        Term newp = (Term) obj;
        return (user == newp.user && m == newp.m);

    }
    @Override
    public String toString(){
        return "(Power :"+String.valueOf(this.p) + ",Rate :"+String.valueOf(this.r)+ ",User :"+String.valueOf(this.user)+ ",M :"+String.valueOf(this.r)+ ",Channel : "+String.valueOf(this.n)+")";
    }
}