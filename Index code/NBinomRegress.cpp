// Model for whales
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
    // Read data
    DATA_MATRIX(X); // Data matrix
    DATA_MATRIX(Xd); // Dispersion matrix
    DATA_MATRIX(Xhat); // Prediction matrix
    DATA_INTEGER(Nyr); // Number of years
    DATA_VECTOR(P_t); // Probability whale remains in area
    DATA_VECTOR(yobs);
    DATA_INTEGER(Eday); // End day for ingress/egress
    DATA_INTEGER(Iday); // Day of ingress
    
    // Parameters
    PARAMETER_VECTOR(beta); // Expected value parameters
    PARAMETER_VECTOR(betad); // dispersion parameters
    
    // Variables
    vector<Type> eta = X*beta; // Expected value parameters
    vector<Type> etad = Xd*betad; // Dispersion
    
    // Run first stage estimation
    Type jnll = 0;
    for (int i=0; i < yobs.size(); i++){
        jnll -= dnbinom_robust(yobs(i), eta(i), 2. * eta(i) - etad(i), true);
        // jnll -= dnbinom2(yobs(i), eta(i), eta(i) * (1 + etad(i)), true);
    }
    
    // Run second stage calculations
    vector<Type> W_t = exp(Xhat*beta); // Estimated number of whales on day t: W_0 = 0
    vector<Type> dW_t = W_t.size(); dW_t.setZero(); // Difference in whales ingressing on day t: W_0 = 0
    vector<Type> A_xyLong = W_t.size(); A_xyLong.setZero(); // Estimated number of whales on day t: W_0 = 0
    vector<Type> A_xy(Nyr); A_xy.setZero();
    Iday = Iday - 1; // Start of ingress - change for cpp indexing
    Eday = Eday - 1; // stop of ingress/egress - change for cpp indexing
    int Nday = 365; // Number of days
    Type tmp_sum = 0;
    
    
    // -- Should be as follows:
    // -- A_0 = 0
    // -- A_1 = dW_1 + p1 * A_0
    // -- A_2 = dW_2 + p1 * A_1 + p2 * A_0
    // -- A_3 = dW_3 + p1 * A_2 + p2 * A_1 + p3 * A_0
    for(int yr = 0; yr < Nyr; yr++){ // Year loop
        for(int k = 0; k < Nday; k++){ // Day loop - Start at day 100 because W_t = 0 for days 0-99. Indexing starts at 0
            tmp_sum = 0; // Initialize
            
            // No whales before April (day 100, index 99)
            if(k < Iday){
                W_t[Nday*yr + k] = 0;
                dW_t[Nday*yr + k] = 0;
                A_xyLong[Nday*yr + k] = 0;
            }
            
            // Accumulation
            if(k >= Iday & k <= Eday){
                dW_t[Nday*yr + k] = W_t[Nday*yr + k] - W_t[Nday*yr + k-1]; // Calculate dW_t
                
                // Run the summation of all previous days (no way to vectorize)
                for(int ktmp = 0; ktmp < k; ktmp++){
                    tmp_sum += dW_t[Nday*yr + ktmp] + P_t(ktmp) * A_xyLong[Nday*yr + k - ktmp - 1];
                }
                
                A_xyLong[Nday*yr + k] = tmp_sum + dW_t[Nday*yr + k]; // segment of k elements, starting at year
            }
            
            // No whales ingressing after day 320 (index 319)
            if(k > Eday){
                dW_t[Nday*yr + k] = 0;
                A_xyLong[Nday*yr + k] = A_xyLong[Nday*yr + Eday];
            }
            
            // Save value
            if(k == Eday){
                A_xy[yr] = A_xyLong[Nday*yr + Eday];
            }
        }
    }  
    
    // Report predictions
    REPORT( W_t );
    REPORT( dW_t );
    REPORT( A_xyLong );
    REPORT( A_xy );
    ADREPORT( A_xy );
    
    return jnll;
}
