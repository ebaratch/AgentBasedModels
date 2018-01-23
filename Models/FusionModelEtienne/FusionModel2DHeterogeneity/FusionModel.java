package Models.FusionModelEtienne.FusionModel2DHeterogeneity;

import Framework.Extensions.SphericalAgent2D;
import Framework.GridsAndAgents.AgentGrid2D;
import Framework.Gui.Vis2DOpenGL;
import Framework.Tools.BitGenome;
import Framework.Utils;

import java.io.File;
import java.io.FileWriter;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Random;

import static Framework.Utils.*;

class Dish extends AgentGrid2D<Cell> {
    final static int BLACK=RGB(0,0,0),RED=RGB(1,0,0),GREEN=RGB(0,1,0),YELLOW=RGB(1,1,0),BLUE=RGB(0,0,1),WHITE=RGB(1,1,1);
    //GLOBAL CONSTANTS
    double DIVISION_PROB=0.0067;
    double DEATH_PROB=0.005;
    double FUSION_PROB=0.05;
    //double FUSION_PROB=0.02;
    //double FUSION_PROB=0.1;

    //BitSet bs = new BitSet();
    //long bits=1+1<<1;
    double CELL_RAD=0.3;
    double MAX_RAD=Math.sqrt(2)*CELL_RAD;
    double FRICTION=0.9;
    double STEADY_STATE_FORCE=0;
    double MAX_STEADY_STATE_LOOPS=10;
    double DIV_RADIUS=CELL_RAD*(2.0/3.0);
    double FORCE_EXPONENT=2;
    double FORCE_SCALER=0.7;
    int GENE_NUMBER=20;
    double mutationProb=0.2;
    int[] WILDTYPE={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int MAXSCORE;
    int[] TypeRepresentants;

    //double MAX_VEL=1000000000;

    int fusionCt=0;
    //INTERNAL VARIABLES
    Random rn=new Random();
    ArrayList<Cell> cellScratch=new ArrayList<>();
    double[] divCoordScratch=new double[2];

    public Dish(int sideLen,int startingPop,double startingRadius){
        super(sideLen,sideLen,Cell.class);
        double[] startCoords=new double[2];
        //for (int i=0;i<GENE_NUMBER;i++) {
        //    WILDTYPE[i]=0;
        //}
        MAXSCORE=0;
        for (int i=0;i<GENE_NUMBER;i++) {
            MAXSCORE=MAXSCORE+(int)Math.pow(2,i);
        }
        TypeRepresentants=new int[MAXSCORE];
        for (int i=0;i<MAXSCORE;i++) {
            TypeRepresentants[i]=0;
        }
        for (int i = 0; i < startingPop; i++) {
            Utils.RandomPointInCircle(startingRadius, startCoords, rn);
            Cell c=NewAgentPT(startCoords[0]+xDim/2.0,startCoords[1]+yDim/2.0);
            if(i%2==0) {
                //for (int j=0;j<GENE_NUMBER;j++) {
                //        WILDTYPE[j]=0;
                //    }
                c.Init(RED,WILDTYPE);

                //c.Mutate();
                //for (int j=0;j<GENE_NUMBER;j++){
                //    System.out.println("WILDTYPE="+j+WILDTYPE[j]);
                //}
            }
            else {
                //for (int j=0;j<GENE_NUMBER;j++) {
                //    WILDTYPE[j]=0;
                //}
                c.Init(RED,WILDTYPE);
            }
        }
    }


    //int[] Blendind(int[] Genotype1,int[] Genotype2){
    //    int[] BlendedGenome=new int[GENE_NUMBER];
    //    for (int i=0;i<GENE_NUMBER;i++){
    //        if(rn.nextDouble()<0.5) {
    //            BlendedGenome[i]=Genotype1[i];
    //            System.out.println("vrai");
    //        }
    //        else{
    //            BlendedGenome[i]=Genotype2[i];
    //            System.out.println("faux");
    //        }
    //    }
    //    return BlendedGenome;
    //}


    void Blendind(int[] Genotype1,int[] Genotype2){
        for (int i=0;i<GENE_NUMBER;i++){
            System.out.println("entries="+i+Genotype1[i]+Genotype2[i]);
        }
        //int[] ScratchGenome=new int[GENE_NUMBER];
        int tempBit=0;
        for (int i=0;i<GENE_NUMBER;i++){
            if(rn.nextDouble()<0.5) {
                tempBit = Genotype1[i];
                //System.out.println("vrai");
                Genotype1[i] = Genotype2[i];
                Genotype2[i] = tempBit;
            }
            else{
                //System.out.println("faux");
            }
        }
        for (int i=0;i<GENE_NUMBER;i++){
            System.out.println("results="+i+Genotype1[i]+Genotype2[i]);
        }
    }


    void Fusion(Cell c1,Cell c2){
        int[] ScratchGenome1=c1.Genotype;
        int[] ScratchGenome2=c2.Genotype;
        //Cell H=NewAgentPT((c1.Xpt()+c2.Xpt())/2,(c1.Ypt()+c2.Ypt())/2);
        //H.Init(c1.color,Blendind(c1.Genotype,c2.Genotype));
        Blendind(ScratchGenome1,ScratchGenome2);

        double X1=c1.Xpt();
        double X2=c2.Xpt();
        double Y1=c1.Ypt();
        double Y2=c2.Ypt();

        c1.Dispose();
        c2.Dispose();

        Cell H1=NewAgentPT(X1,Y1);
        Cell H2=NewAgentPT(X2,Y2);

        H1.Init(RED,ScratchGenome1);
        H2.Init(RED,ScratchGenome2);

    }

    int SteadyStateMovement(){
        int loopCt=0;
        while(loopCt<MAX_STEADY_STATE_LOOPS) {
            double maxForce=0;
            for (Cell c : this) {
                maxForce=Math.max(maxForce,c.Observe());
            }
            for (Cell c : this) {
                c.Act();
            }
            loopCt++;
            if(maxForce<STEADY_STATE_FORCE){
                //System.out.println(loopCt+","+maxForce);
                break;
            }
        }
        return loopCt;
    }
    void Step(){
        SteadyStateMovement();
        //for (Cell c:this) {
        //    c.Observe();
        //}
        //for (Cell c:this){
        //    c.Act();
        //}
        for (Cell c:this) {
            c.Step();
        }
        for (Cell c:this){
            fusionCt+=c.Fuse()?1:0;
        }
        //System.out.println(fusionCt);
        IncTick();
    }


    int ConvertBinary(int[] BinaryVector){
        int total=0;
        for (int i=0;i<BinaryVector.length;i++){
            total=total+BinaryVector[i]*(int)Math.pow(2,i);
        }
        return total;
    }


    double ComputeShannon(){
        double Shanon=0;
        int GenoInt=0;
        int TotalCell=0;
        double sumProp=0;

        for (Cell c:this){
            TotalCell=TotalCell+1;
            GenoInt=ConvertBinary(c.Genotype);
            TypeRepresentants[GenoInt]=TypeRepresentants[GenoInt]+1;
        }
        for (int i=0;i<MAXSCORE;i++){
            if(TypeRepresentants[i]>0){
                Shanon=Shanon-((float)TypeRepresentants[i]/TotalCell)*Math.log((float)TypeRepresentants[i]/TotalCell);
                System.out.println("checkProp="+i+(float)TypeRepresentants[i]/TotalCell);
                sumProp=sumProp+(float)TypeRepresentants[i]/TotalCell;
            }
            TypeRepresentants[i]=0;
        }
        System.out.println("checkPop="+TotalCell+","+sumProp);
        return Shanon;
    }
}

class Cell extends SphericalAgent2D<Cell,Dish> {
    int color;
    boolean hybrid;
    int[] Genotype;


    void Init(int InitialColor,int[] GenotypeIni){
        radius=G().CELL_RAD;
        xVel=0;
        yVel=0;
        //color=InitialColor;
        hybrid=false;
        Genotype=new int[G().GENE_NUMBER];
        for (int i=0;i<G().GENE_NUMBER;i++){
            Genotype[i]=GenotypeIni[i];
        }
        //for (int i=0;i<G().GENE_NUMBER;i++){
        //    Genotype[i]=0;
        //}
        SetCellColor();
    }
    void Init(int InitialColor,boolean IsHybrid,int[] GenotypeIni){
        xVel=0;
        yVel=0;

        //Genotype=new int[G().GENE_NUMBER];
        Genotype=GenotypeIni;
        //for (int i=0;i<G().GENE_NUMBER;i++){
        //    Genotype[i]=0;
        //}
        hybrid=IsHybrid;
        if(hybrid==false){
            radius=G().CELL_RAD;
        }
        else{
            radius=Math.sqrt(2)*G().CELL_RAD;
        }
        //color=InitialColor;
        SetCellColor();
    }

    //void SetCellColor(int newColor){
    //    color=newColor;
    //}

    void SetCellColor(){
        double sum=0;
        for (int i=0;i<G().GENE_NUMBER;i++){
            sum=sum+Genotype[i];
        }
        color = HeatMapRGB(((sum)+2)/G().GENE_NUMBER);
        //if (sum==0){
        //    color=G().GREEN;
        //}
        //else if (sum==1.0){
        //    color=G().RED;
        //}
        //else if (sum==2.0){
        //    color=G().BLUE;
        //}
    }


    void Mutate(){
        for (int i=0;i<G().GENE_NUMBER;i++){
            System.out.println("entries="+i+Genotype[i]);
        }
        int indexMut=G().rn.nextInt(G().GENE_NUMBER);
        if (Genotype[indexMut]==0){
            Genotype[indexMut]=1;
        }
        else{
            Genotype[indexMut]=0;
        }
        for (int i=0;i<G().GENE_NUMBER;i++){
            System.out.println("sorties="+i+Genotype[i]);
        }
        SetCellColor();
    }

    double OverlapToForce(double overlap){
        if(overlap<0){
            return 0;
        }
        return Math.pow(G().FORCE_SCALER*overlap,G().FORCE_EXPONENT);
        //return G().FORCE_SCALER*overlap;
    }
    boolean Fuse(){
        if(hybrid){return false;}
        //listing all cells in the area
        G().cellScratch.clear();
        G().AgentsInRad(G().cellScratch,Xpt(),Ypt(),G().CELL_RAD*2);
        int neighborCt=0;
        //getting valid fusion neighbors
        for (int i=0;i<G().cellScratch.size();i++) {
            Cell c=G().cellScratch.get(i);
            if(!c.hybrid&&c!=this&&Utils.DistSquared(Xpt(),Ypt(),c.Xpt(),c.Ypt())<G().CELL_RAD*2){
                G().cellScratch.set(neighborCt,c);
                neighborCt++;
            }
        }
        //fusing
        if(neighborCt>0&&G().rn.nextDouble()<Utils.ProbScale(G().FUSION_PROB,neighborCt)){
            G().Fusion(this,G().cellScratch.get(G().rn.nextInt(neighborCt)));
            return true;
        }
        return false;
    }
    double Observe(){

        return SumForces(radius+G().MAX_RAD,G().cellScratch,this::OverlapToForce);
        /*
        if(ret>G().MAX_VEL){
            xVel*=G().MAX_VEL/ret;
            yVel*=G().MAX_VEL/ret;
        }
        return ret;
        */
    }
    void Act(){
        ForceMove();
        ApplyFriction(G().FRICTION);
    }
    void Step(){
        if(G().rn.nextDouble()<G().DEATH_PROB && hybrid==false){
            Dispose();
            return;
        }
        if(G().rn.nextDouble()<G().DIVISION_PROB && hybrid==false){
            Cell child=Divide(G().DIV_RADIUS,G().divCoordScratch,G().rn);
            if(G().rn.nextDouble()<G().mutationProb){
                child.Init(this.color,this.Genotype);
                child.Mutate();
                //child.MutColorChange(this,G().rn,0.1);
            }
            else{
                child.Init(this.color,this.Genotype);
            }
        }

    }
}

public class FusionModel {
    static int SIDE_LEN = 100;
    static int STARTING_POP = 1000;
    static double STARTING_RADIUS = 10;
    static int TIMESTEPS = 2000;
    static float[] circleCoords = Utils.GenCirclePoints(1, 10);

    public static void main(String[] args) {
        //TickTimer trt=new TickRateTimer();

        Vis2DOpenGL vis = new Vis2DOpenGL("Cell Fusion Visualization", 1000, 1000, SIDE_LEN, SIDE_LEN);
        Dish d = new Dish(SIDE_LEN, STARTING_POP, STARTING_RADIUS);
        double Shanon=0;
        //double Shanon = 0;
        //int[] vectorExample={1,0,0,0,1,0};
        //int[] vectorExample2={0,1,1,1,1,0};
        //for (int i=0;i<6;i++){
        //    System.out.println("entries="+i+vectorExample[i]+vectorExample2[i]);
        //}
        //d.Blendind(vectorExample,vectorExample2);
        //for (int i=0;i<6;i++){
        //    System.out.println("results="+i+vectorExample[i]+vectorExample2[i]);
        //}

        //d.SetCellsColor("red");



        //for (int i = 0; i < TIMESTEPS; i++) {
        //    vis.TickPause(0);
        //    d.Step();
        //    DrawCells(vis,d,"/Users/baratcEA/work/Moffitt/Phoenix/Fusion/ResultsSpatial/Heterogeneity2D/Mutation/",i);
        //    Shanon=0;
        //    Shanon=d.ComputeShannon();

        //    System.out.println("Shanon index=" + Shanon);
        //}

        try {
            File file = new File("/Users/baratcEA/work/Moffitt/Phoenix/Fusion/ResultsSpatial/Heterogeneity2D/Mutation/Shanon.txt");
            PrintWriter printWriter = new PrintWriter(file);

            for (int i = 0; i < TIMESTEPS; i++) {
                vis.TickPause(0);
                d.Step();

                DrawCells(vis, d, "/Users/baratcEA/work/Moffitt/Phoenix/Fusion/ResultsSpatial/Heterogeneity2D/Mutation/", i);

                Shanon = d.ComputeShannon();

                String line = "";
                line = "" + i + "," + Shanon + "\n";
                System.out.println("Shanon index=" + Shanon);
                printWriter.print(line);
            }

            printWriter.close();
        }// end try block
        catch (Exception e) {
            System.out.println(e.getClass());
        }
    }



    static void DrawCells(Vis2DOpenGL vis,Dish d,String path,int i){
        vis.Clear(Dish.WHITE);
        for (Cell c:d) {
            //color "cytoplasm"
            vis.Circle(c.Xpt(),c.Ypt(),c.radius,c.color);
        }
        for (Cell c:d) {
            //color "nucleus"
            vis.Circle(c.Xpt(),c.Ypt(),c.radius/2.0, Dish.BLUE);
        }
        vis.Show();
        //vis.ToPNG((path.concat(Integer.toString(i))).concat(".png"));
    }
}
