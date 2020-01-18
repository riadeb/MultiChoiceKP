import java.util.ArrayList;

public class pair {
    int p;
    int r;
    int user; // user index
    int m; //power index
    int n; //channel index
    double inc_eff; //equals inc_rate/inc_eff
    int inc_power;
    int inc_rate;
    public  pair(int p, int r, int u, int m, int i) {
        this.p = p;
        this.r= r;
        this.user = u;
        this.m = m;
        this.n = i;
    }

    @Override
    public boolean equals(Object obj) {
        pair newp = (pair) obj;
        return (user == newp.user && m == newp.m);

    }
    @Override
    public String toString(){
        return "(Power :"+String.valueOf(this.p) + ",Rate :"+String.valueOf(this.r)+ ",User :"+String.valueOf(this.user)+ ",M :"+String.valueOf(this.r)+ ",Channel : "+String.valueOf(this.n)+")";
    }
}

class sol_pair extends pair {
    double x;
    public sol_pair(int p, int r, int u, int m,int i, double x) {
        super(p,r,u,m,i);
        this.x = x;
    }
}

class Solution {
    double Rate;
    ArrayList<sol_pair>[] data;
    public Solution(double r, ArrayList<sol_pair>[] d) {
        Rate = r;
        data = d;
    }
}