import java.lang.reflect.Array;
import java.util.*;

public class onlineAlg{
    int p;
    int pmax;
    int rmax;
    int M;
    int N;
    int K;
    

    public onlineAlg(int p,int pmax,int rmax,int M,int N,int K)
    {
        this.p = p;
        this.pmax = pmax;
        this.rmax = rmax;
        this.M = M;
        this.N = N;
        this.K = K;
    }

    public static void main(String[] args) throws Exception {
       
        onlineAlg online = new onlineAlg(100, 50, 100, 2, 4, 10);
        float mean = 0;
        float meanp = 0;
        float meanpdp = 0;
        for (int i = 0;i<10000;i++){
            dt Dt = online.generate();
            int[] e = online.StochasticOnlineScheduling(Dt.data);
            float a = e[0];
            meanp += e[1];
            int[] f = online.DP(Dt.Data);
            meanpdp += f[1];
            float b = f[0];
            mean += a/b;
        }
        System.out.println(mean/10000);
        System.out.println(meanp/10000);
        System.out.println(meanpdp/10000);
        


    }
    dt generate(){
        ArrayList<pair>[] data = new ArrayList[this.K];
        ArrayList<pair>[] Data = new ArrayList[this.N];

        for (int k = 0; k < this.K; k++) {
            data[k] = new ArrayList<pair>();
        }
        for (int n = 0; n < N; n++) {
            Data[n] = new ArrayList<pair>();
        }
        
        Random rd = new Random();
        
        for (int k = 1; k <= this.K; k++) {
            for (int n = 1; n <= this.N; n++) {
                for (int m = 1; m <= this.M; m++) {
                    int pRandom = rd.nextInt(this.pmax) + 1;
                    int rRandom = rd.nextInt(this.rmax) + 1;
                    
                    pair pr = new pair(pRandom, rRandom, k, m , n );
                    data[k - 1].add(pr);
                    Data[n-1].add(pr);
                }
            }
        }
        dt Dt = new dt(data,Data);
        
        return Dt;
    }
    int[] StochasticOnlineScheduling(ArrayList<pair>[] data) {
        
        int powerBud = this.p;
        int achievedRate = 0;
        double expectation = 0;
        int[] listOfChannels = new int[this.N];
        Random rd = new Random();
        int ChannelPerUser = (int) this.N / this.K;
        int PowerPerUser = (int) this.p / this.K;
        for (int j = 0;j<this.N;j++){
            listOfChannels[j] = 0;
        }
        for (int k = 1; k <= this.K; k++) {
            
            int countChannels = 0;
            
            
            Collections.sort(data[k-1], new eff_comparator());
            int up = 0;
            int i =0;
            while (up < PowerPerUser && i<data[0].size()) {
                pair pai = data[k-1].get(i);
                if (listOfChannels[pai.n-1] == 0 && up+pai.p <=PowerPerUser){
                    listOfChannels[pai.n-1] =1;
                    countChannels +=1;
                    powerBud-= pai.p;
                    up+=pai.p;
                    achievedRate+=pai.r;
                }
                
                i+=1;
                
            }
            up = 0;
            

        }
        
        int[] ans = new int[2];
        ans[0]= achievedRate;
        ans[1] = powerBud;
        return ans;
    }

    int[] DP(ArrayList<pair>[] Data){ //First Implementation of Dynamic programming
        /*
        L[p-1] holds maximum rate with power budget p
          */
        //System.out.println(Data[0]);
        int[] L = new int[p];
        for(int i = 1; i <= p; i++) {
            int maxr = 0;
            for(pair curr_pair:Data[0]){
                if (curr_pair.p <= i) maxr = Math.max(maxr,curr_pair.r);
            }
            L[i-1] = maxr;
            
        }
        for(int i =1; i < N; i++){
            int[] tempL = new int[p];
            for (int power = 1; power <= p; power++) {
                ArrayList<pair> curr_channel = Data[i];
                int max_r = 0;
                for(pair pp:curr_channel) {
                    if (pp.p < power && L[power-pp.p - 1] > 0) max_r = Math.max(max_r,L[power-pp.p - 1] + pp.r);
                }
                tempL[power-1] = max_r;
            }
            L = tempL;
            //System.out.print("######n");System.out.print(i);System.out.println(", ");
            /*for (int lk:L){
                System.out.print(lk);System.out.print(", ");
            }*/
        }
        int b = 0;
        int po = L[p-1];
        while (L[p-2-b] == po ){
            b+=1;
        }
        int[] ans = new int[2];
        ans[0] = L[p-1];
        ans[1] = p-b;
        return ans;
    }

}
class eff_comparator implements  Comparator {
    public int compare(Object o1, Object o2) {
        pair pair1 = (pair) o1;
        pair pair2 = (pair) o2;
        return Double.compare(pair1.r/pair1.p,pair2.r/pair2.p);
    }
}
class dt {
    ArrayList<pair>[] data;
    ArrayList<pair>[] Data;
    public dt(ArrayList<pair>[] data, ArrayList<pair>[] Data){
        this.data = data;
        this.Data = Data;
    }
}