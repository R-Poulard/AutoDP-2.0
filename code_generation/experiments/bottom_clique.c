
int C(int i,int j,int k, int l){
    return 0;
}
int C2(int i,int j, int k, int l){
    return 0;
}

//important aspect is the top to bottom situation (calculating from tree leaf to root)

int main(){

    //CLIQUE CASE
    int ianch_max=10;
    int ianch_min=5;


    int janch_max=10;
    int janch_min=5;
    int j=0;
    int k=0;
        // i goes in reverse to get i+1 values
        //j goes forward to get j-1 values
        // anch should be i l

        //maybe good to put the get inside the function (no call back to C(j-1 ect))
        //but then the INT_MAX 'sanity check' has to be here (otherwise get value who doesn't exist)
    for(int i=k;i>ianch_min;i--){
        for(int l=janch_max;l<=j;l++){
            C(i,j,k,l);
            C2(i,j,k,l);
        }
    }

    //Transitional bag is already top to bottom i feel
        //but changes to "get" inside the loop not calling a function
    
    //Diagonal    
        // sum in D should also be replaced by get as it as already been calculated
        // maybe not the best as sanity check would also has to be performed here (add modularity)

    //here the boudaries would be S values surrounding i and j
    //if get inside sanity check needs to be done here also
    for(int i=k;i>ianch_min;i--){
        for(int j=janch_max;j<=j;l++){
            D(i,j,l);
            D2(i,j,l);
        }
    }
}

//fold=>
//calculate_Z()
//calculate_X()
//calculate_Y()
//...
//score=calculate_A()

// backtrack?
// same thing exactly it feels?

