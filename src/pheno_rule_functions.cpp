# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List inflection_cpp(arma::vec y){
  if (y.n_elem < 3) {
    stop("y must contain more than 2 elements.");
  }
  arma::vec dif_y = y.subvec(1,y.n_elem - 1) - y.subvec(0,y.n_elem - 2);
  arma::vec s_y(dif_y.n_elem);
  s_y.zeros();
  arma::uvec pos = find(dif_y > 0);
  s_y(pos).fill(1.0);
  s_y(find(dif_y < 0)).fill(-1.0);
  arma::vec id = s_y.subvec(1,s_y.n_elem - 1) - s_y.subvec(0,s_y.n_elem - 2);
  arma::vec time = conv_to<arma::vec>::from(find(id != 0)) + 2;
  arma::vec type(time.n_elem);
  type.ones();
  arma::vec sos = conv_to<arma::vec>::from(find(id > 0)) + 2;
  for(int i=0; i < type.n_elem; i++){
    if(any(sos == time(i))){
      type(i) = 2;
    }
  }
  arma::mat df(type.n_elem, 2);
  df.col(0) = time;
  df.col(1) = y(find(id!=0)+1);
  
  return List::create(
    Named("df") = df,
    Named("type") = type);
}

// [[Rcpp::export]]
int armaWhichMin(arma::vec x){
  return which_min(as<NumericVector>(wrap(x)));
}

// [[Rcpp::export]]
int armaWhichMax(arma::vec x){
  return which_max(as<NumericVector>(wrap(x)));
}

// [[Rcpp::export]]
List remove_sos_new(arma::mat pog_dat, arma::mat sos_dat) {
  // pog_dat must have at least one row
  // sos_dat must have at least one row
  arma::uvec sos_keep(sos_dat.n_rows);
  sos_keep.zeros();
  arma::vec tpk = pog_dat.col(0);
  arma::vec sos_time = sos_dat.col(0);
  arma::vec sos_level = sos_dat.col(1);
  
  int p = tpk.n_elem;
  arma::uvec s = find(sos_time >= 1 && sos_time < tpk(0));
  
  if(s.n_elem != 0){
    sos_keep(s(0) + armaWhichMin(sos_level.elem(s))) = 1;
  }
  for(int i=0; i < (p - 1); i++){
    s = find(sos_time > tpk(i) && sos_time < tpk(i+1));
    if(s.n_elem != 0){
      sos_keep(s(0) + armaWhichMin(sos_level.elem(s))) = 1;
    }
  }
  s = find(sos_time > tpk(p-1) && sos_time <= 230);
  if(s.n_elem != 0){
    sos_keep(s(0) + armaWhichMin(sos_level.elem(s))) = 1;
  }
  
  arma::mat ss = sos_dat.rows(find(sos_keep==1));
  
  return List::create(
    Named("pog_dat") = pog_dat,
    Named("sos_dat") = ss);
}

// [[Rcpp::export]]
List remove_pog_new(arma::mat sos_dat, arma::mat pog_dat) {
  // sos_dat must have at least one row
  // pog_dat must have at least one row

  arma::uvec pog_keep(pog_dat.n_rows);
  pog_keep.zeros();
  arma::vec pog_time = pog_dat.col(0);
  arma::vec tsk = sos_dat.col(0);
  arma::vec pog_level = pog_dat.col(1);

  // keep largest POG between kept SOS
  int p = tsk.n_elem;
  arma::uvec s = find(pog_time >= 1 && pog_time < tsk(0));
  if(s.n_elem != 0){
    pog_keep(s(0) + armaWhichMax(pog_level.elem(s))) = 1;
  }
  for(int i=0; i < (p - 1); i++){
    s = find(pog_time > tsk(i) && pog_time < tsk(i+1));
    if(s.n_elem != 0){
      pog_keep(s(0) + armaWhichMax(pog_level.elem(s))) = 1;
    }
  }

  s = find(pog_time > tsk(p-1) && pog_time <= 230);
  if(s.n_elem != 0){
    pog_keep(s(0) + armaWhichMax(pog_level.elem(s))) = 1;
  }

  arma::mat pp = pog_dat.rows(find(pog_keep==1));

  return List::create(
    Named("pog_dat") = pp,
    Named("sos_dat") = sos_dat);
}

// [[Rcpp::export]]
List merge_df2(arma::mat sos_dat, arma::mat pog_dat) {
  int sl = sos_dat.n_rows;
  int pl = pog_dat.n_rows;
  arma::mat df(sl+pl, 2);
  NumericVector ss, ps;
  if(sos_dat(0,0) < pog_dat(0,0)){
    ss = seq(0,sl-1)*2.0;
    ps = seq(0,pl-1)*2.0+1;
  }else{
    ss = seq(0,sl-1)*2.0+1;
    ps = seq(0,pl-1)*2.0;
  }
  NumericVector type(sl + pl);
  NumericVector pog(pl, 1.0);
  NumericVector sos(sl, 2.0);

  type[ps] = pog;
  type[ss] = sos;
  for(int i = 0; i < pl; i++){
    df.row(ps(i)) = pog_dat.row(i);
  }
  for(int i = 0; i < sl; i++){
    df.row(ss(i)) = sos_dat.row(i);
  }
  
  if(df.n_rows < 2){
    return List::create(
      Named("df") = arma::mat(1,2).fill(NA_REAL),
      Named("type") = NA_REAL);
  }else{
    return List::create(
      Named("df") = df,
      Named("type") = type);
  }
}

// [[Rcpp::export]]
List checkd_min_new(arma::mat dat, arma::vec type, int d) {
  arma::mat df_pog = dat.rows(find(type==1));
  arma::mat df_sos = dat.rows(find(type==2));
  
  if(df_pog.n_rows < 2){
    stop("Need at least two POG for checkd_min_new()");
  }
  
  arma::uvec keep(df_pog.n_rows);
  keep.ones();
  
  arma::vec pog_time = df_pog.col(0);
  arma::vec dif_time = pog_time.subvec(1,df_pog.n_rows - 1) - pog_time.subvec(0,df_pog.n_rows - 2);
  
  int ds = dif_time.n_elem;
  arma::vec level = df_pog.col(1);
  
  for (int i=0; i < ds; i++) {
    if(dif_time(i) < d){
      if(level(i) < level(i+1)){
        keep(i) = 0;
      }else{
        keep(i+1) = 0;
      }
    }
  }
  
  arma::mat df_pog2 = df_pog.rows(find(keep==1));
  
  if(df_pog2.n_rows >= 2){
    List t2 = remove_sos_new(df_pog2, df_sos);
    return merge_df2(t2[1], t2[0]);
  }else if(df_pog2.n_rows==1 && df_pog2(0,0) > df_sos(0,0)){
    List t2 = remove_pog_new(df_sos, df_pog2);
    return merge_df2(t2[1], t2[0]);
  }else{
    return List::create(
      Named("df") = arma::mat(1,2).fill(NA_REAL),
      Named("type") = NA_REAL);
  }
}

// [[Rcpp::export]]
List checkd_max_new(arma::mat dat, arma::vec type, int d) {
  arma::mat df_pog = dat.rows(find(type==1));
  arma::mat df_sos = dat.rows(find(type==2));
  
  if(df_sos.n_rows < 2){
    stop("Need at least two SOS for checkd_max_new()");
  }
  
  arma::uvec keep(df_sos.n_rows);
  keep.ones();
  
  arma::vec sos_time = df_sos.col(0);
  arma::vec dif_time = sos_time.subvec(1,df_sos.n_rows - 1) - sos_time.subvec(0,df_sos.n_rows - 2);
  
  int ds = dif_time.n_elem;
  arma::vec level = df_sos.col(1);
  
  for (int i=0; i < ds; i++) {
    if(dif_time(i) < d){
      if(level(i) < level(i+1)){
        keep(i+1) = 0;
      }else{
        keep(i) = 0;
      }
    }
  }
  
  arma::mat df_sos2 = df_sos.rows(find(keep==1));
  
  if(df_sos2.n_rows >= 2){
    List t2 = remove_pog_new(df_sos2, df_pog);
    return merge_df2(t2[1], t2[0]);
  }else if(df_sos2.n_rows==1 && df_pog(0,0) > df_sos2(0,0)){
    List t2 = remove_pog_new(df_sos2, df_pog);
    return merge_df2(t2[1], t2[0]);
  }else{
    return List::create(
      Named("df") = arma::mat(1,2).fill(NA_REAL),
      Named("type") = NA_REAL);
  }
}

// [[Rcpp::export]]
List checkd_max_alt(arma::mat dat, arma::vec type, int d) {
  arma::mat df_pog = dat.rows(find(type==1));
  arma::mat df_sos = dat.rows(find(type==2));
  
  if(df_sos.n_rows < 2){
    stop("Need at least two SOS for checkd_max_alt()");
  }
  if(df_pog.n_rows < 1){
    stop("Need at least one POG for checkd_max_alt()");
  }
  
  arma::uvec keep(df_sos.n_rows);
  keep.ones();
  
  int ds = df_sos.n_rows;
  arma::vec times = df_sos.col(0);
  arma::vec level = df_sos.col(1);
  int i = 0, j = 1;
  while (i < ds && j < ds) {
    if(times(j) - times(i) < d){
      if(level(i) < level(j)){
        keep(j) = 0;
        j++;
      }else{
        keep(i) = 0;
        i=j;
        j++;
      }
    }else{
      i=j;
      j++;
    }
  }
  
  arma::mat df_sos2 = df_sos.rows(find(keep==1));
  
  if(df_sos2.n_rows >= 2){
    List t2 = remove_pog_new(df_sos2, df_pog);
    return merge_df2(t2[1], t2[0]);
  }else if(df_sos2.n_rows==1 && df_pog(0,0) > df_sos2(0,0)){
    List t2 = remove_pog_new(df_sos2, df_pog);
    return merge_df2(t2[1], t2[0]);
  }else{
    return List::create(
      Named("df") = arma::mat(1,2).fill(NA_REAL),
      Named("type") = NA_REAL);
  }
}

// [[Rcpp::export]]
List check_2sos(arma::mat df, arma::vec type, double cf){
  // function to check sos when diff between previous and next pog are both < cf
  // must have at least 1 element
  if(type.n_elem<=1){
    return List::create(
      Named("df") = arma::mat(1,2).fill(NA_REAL),
      Named("type") = NA_REAL);
  }else{
    arma::uvec keep(df.n_rows);
    keep.ones();
    arma::uvec sos = find(type == 2);
    arma::mat df_sos, df_pog;
    int n = df.n_rows;
    for(int i=0; i < sos.n_elem; i++){ // for each of i SOS
      if(sos(i) > 1 && sos(i) < n-1){ // must be at least two event prior and two events past the considered SOS
        if(std::abs(df(sos(i)+1,1) - df(sos(i),1)) < cf && std::abs(df(sos(i)-1,1) - df(sos(i),1)) < cf){
          // check if this difference in previous pog and sos(i) < cf and next pog and sos(i) < cf
          if(df(sos(i),1) > df(sos(i)-2,1)){
            // then if current sos(i) > last sos(i), remove current sos(i)
            keep(sos(i)) = 0;
          }
          if(df(sos(i)-1,1) < df(sos(i)+1,1)){
            // then check which of surrounding pog is larger, remove smaller one
            keep(sos(i)-1) = 0;
          }else{
            keep(sos(i)+1) = 0;
          }
        }
      }
    }
    arma::mat df2 = df.rows(find(keep==1));
    arma::vec type2 = type(find(keep==1));
    df_sos = df2.rows(find(type2==2));
    df_pog = df2.rows(find(type2==1));
    
    if(df_sos.n_rows < 1 || df_pog.n_rows < 1){
      return List::create(
        Named("df") = arma::mat(1,2).fill(NA_REAL),
        Named("type") = NA_REAL);
    }else{
      List t3 = remove_sos_new(df_pog, df_sos);
      df_sos = as<arma::mat>(t3[1]);
      df_pog = as<arma::mat>(t3[0]);
      if(df_sos.n_rows < 1 || df_pog.n_rows < 1){
        return List::create(
          Named("df") = arma::mat(1,2).fill(NA_REAL),
          Named("type") = NA_REAL);
      }else{
        List t4 = remove_pog_new(df_sos, df_pog);
        df_sos = as<arma::mat>(t4[1]);
        df_pog = as<arma::mat>(t4[0]);
        if(df_sos.n_rows < 1 || df_pog.n_rows < 1){
          return List::create(
            Named("df") = arma::mat(1,2).fill(NA_REAL),
            Named("type") = NA_REAL);
        }else{
          return merge_df2(df_sos, df_pog);
        }
      }
    }
  }
}

// [[Rcpp::export]]
List check_2pog(arma::mat df, arma::vec type, double cf){
  // function to check pog when diff between previous and next sos are both < cf
  // must have at least 1 element
  if(type.n_elem <= 1){
    return List::create(
      Named("df") = arma::mat(1,2).fill(NA_REAL),
      Named("type") = NA_REAL);
  }else{
    arma::uvec keep(df.n_rows);
    keep.ones();
    arma::uvec sos = find(type == 2);
    arma::mat df_sos, df_pog;
    int n = df.n_rows;
    for(int i=0; i < sos.n_elem; i++){
      if(sos(i) < n-2){ // only consider sos that aren't the last two 
        if(sos(i)==0){ // if sos(i) is the first identified event
          if(std::abs(df(sos(i)+1,1) - df(sos(i),1)) < cf && std::abs(df(sos(i)+1,1) - df(sos(i)+2,1)) < cf){
            // sos(i) + 1 is first identified pog, if difference between sos, pog, and sos + 1 < cf
            keep(sos(i)+1) = 0;
          }
        }else if(std::abs(df(sos(i)+1,1) - df(sos(i),1)) < cf || std::abs(df(sos(i)-1,1) - df(sos(i),1)) < cf){
          // current pog and prev sos < cf OR previous pog and previous sos < cf
          if(std::abs(df(sos(i)+1,1) - df(sos(i),1)) < cf && std::abs(df(sos(i)+1,1) - df(sos(i)+2,1)) < cf){
            keep(sos(i)+1) = 0;
            if(df(sos(i),1) < df(sos(i)+2,1)){
              keep(sos(i)+2) = 0;
            }else{
              keep(sos(i)) = 0;
            }
          }
        }
      }
    }
  
    arma::mat df2 = df.rows(find(keep==1));
    arma::vec type2 = type(find(keep==1));
    df_sos = df2.rows(find(type2==2));
    df_pog = df2.rows(find(type2==1));
    
    if(df_sos.n_rows < 1 || df_pog.n_rows < 1){
      return List::create(
        Named("df") = arma::mat(1,2).fill(NA_REAL),
        Named("type") = NA_REAL);
    }else{
      List t3 = remove_sos_new(df_pog, df_sos);
      df_sos = as<arma::mat>(t3[1]);
      df_pog = as<arma::mat>(t3[0]);
      if(df_sos.n_rows < 1 || df_pog.n_rows < 1){
        return List::create(
          Named("df") = arma::mat(1,2).fill(NA_REAL),
          Named("type") = NA_REAL);
      }else{
        List t4 = remove_pog_new(df_sos, df_pog);
        df_sos = as<arma::mat>(t4[1]);
        df_pog = as<arma::mat>(t4[0]);
        if(df_sos.n_rows < 1 || df_pog.n_rows < 1){
          return List::create(
            Named("df") = arma::mat(1,2).fill(NA_REAL),
            Named("type") = NA_REAL);
        }else{
          return merge_df2(df_sos, df_pog);
        }
      }
    }
  }
}

// [[Rcpp::export]]
List check_2up(arma::mat df, arma::vec type, arma::vec y, double cf){
  if(type.n_elem <= 1){
    return List::create(
      Named("df") = arma::mat(1,2).fill(NA_REAL),
      Named("type") = NA_REAL);
  }else{
    arma::uvec keep(df.n_rows);
    keep.ones();
    arma::uvec sos = find(type == 2);
    arma::mat df_sos, df_pog;
    int n = df.n_rows;
    for(int i=0; i < sos.n_elem; i++){
      // if n = 2, never enter if statements
      if(sos(i) < n-2){ // consider sos not the last two events
        if(sos(i)==1){ // if there is a pog before first sos
          if(df(sos(i)-1,1) - df(sos(i),1) < cf && df(sos(i)-1,1) - df(sos(i),1) > 0){
            // if difference between sos and previous pog is < cf but greater than 0
            if(df(sos(i)+1,1) - df(sos(i),1) > cf){
              // if difference between sos and next pog is > cf, remove both sos and prev pog
              keep(sos(i)-1)=0;
              keep(sos(i))=0;
            }
          }
        }else if(sos(i) > 1){
          if(df(sos(i)-1,1) - df(sos(i),1) < cf && df(sos(i)-1,1) - df(sos(i),1) > 0){
            // if prev pog and sos < cf and greater than 0
            if(df(sos(i)-1,1) - df(sos(i)-2,1) > cf && df(sos(i)+1,1) - df(sos(i),1) > cf){
              keep(sos(i)-1)=0;
              keep(sos(i))=0;
              // if difference in previous pog and previous sos > cf, and next pog and current sos > cf, remove current sos and prev pog
            }
          }else if(df(sos(i)+1,1) - df(sos(i),1) < cf && df(sos(i)+1,1) - df(sos(i),1) > 0){
            // if current sos and pog < cf and greater than 0
            if(i==(sos.n_elem-1) && df(sos(i)-1,1) - df(sos(i),1) > cf){
              keep(sos(i)+1)=0;
              keep(sos(i))=0;
              // if sos is second to last and difference between prev pos and sos > cf, get rid of sos and next pog
            }else if(df(sos(i)-1,1) - df(sos(i),1) > cf & df(sos(i)+1,1) - df(sos(i)+2,1) > cf){
              keep(sos(i)+1)=0;
              keep(sos(i))=0;
              // if prev pog and current sos > cf and current pog and next sos > cf, remove current sos and current pog
            }
          }
        }
      }
    }
  
    arma::uvec kp = find(keep==1);
    arma::mat df2 = df.rows(kp);
    arma::vec type2 = type(kp);
    arma::uvec keep2 = keep(kp);
  
    if(df2.n_rows > 1){
    // fix weird last stuff
      int fp = df2(0,0);
    
      if(fp == 1.0){
        keep2(0) = 0;
      }else if(type2(0)==1){
        if(any(y.subvec(0,fp-2) > df2(0,1))){
          keep2(0) = 0;
        }
      }else if(type2(0)==2){
        if(any(y.subvec(0,fp-2) < df2(0,1))){
          keep2(0) = 0;
        }
      }

      if(std::abs(df2(0,1) - df2(1,1)) < cf){
        keep2(0) = 0;
      }
      
      int r = df2.n_rows - 1;
      int lp = df2(r,0);
      int ny = y.n_elem;
      if(lp == ny){
          keep2(r) = 0;
      }else if(type2(r)==1){
        if(any(y.subvec(lp,ny-1) > df2(r,1))){
          keep2(r) = 0;
        }
      }else if(type2(r)==2){
        if(any(y.subvec(lp,ny-1) < df2(r,1))){
          keep2(r) = 0;
        }
      }
    
      if(std::abs(df2(r,1) - df2(r-1,1)) < cf){
        keep2(r) = 0;
      }
    }
    arma::uvec kp2 = find(keep2==1);
    arma::mat df3 = df2.rows(kp2);
    arma::vec type3 = type2(kp2);
  
    if(df2.n_rows <=1 || df3.n_rows <= 1){
      return List::create(
        Named("df") = arma::mat(1,2).fill(NA_REAL),
        Named("type") = NA_REAL);
    }else{
      df_sos = df3.rows(find(type3==2));
      df_pog = df3.rows(find(type3==1));
      if(df_sos.n_rows < 1 || df_pog.n_rows < 1){
        return List::create(
          Named("df") = arma::mat(1,2).fill(NA_REAL),
          Named("type") = NA_REAL);
      }else{
        List df4 = remove_sos_new(df_pog, df_sos);
        df_sos = as<arma::mat>(df4[1]);
        df_pog = as<arma::mat>(df4[0]);
        if(df_sos.n_rows < 1 || df_pog.n_rows < 1){
          return List::create(
            Named("df") = arma::mat(1,2).fill(NA_REAL),
            Named("type") = NA_REAL);
        }else{
          List df5 = remove_pog_new(df_sos, df_pog);
          df_sos = as<arma::mat>(df5[1]);
          df_pog = as<arma::mat>(df5[0]);
          if(df_sos.n_rows < 1 || df_pog.n_rows < 1){
            return List::create(
              Named("df") = arma::mat(1,2).fill(NA_REAL),
              Named("type") = NA_REAL);
          }else{
            return merge_df2(df_sos, df_pog);
          }
        }
      }
    }
  }
}

// [[Rcpp::export]]
List rule_cpp_no_errors(arma::vec y, double cf, int d, int s){
  if(d > s){
    stop("L cannot be greater than the period, s.");
  }
  List inf = inflection_cpp(y); // find all maxima and minima
  arma::mat df = inf[0]; // first column is time, second column is value
  arma::vec type = inf[1]; // vector defining whether maxima/pog (1) or minima/sos (2)
  

  if(df.n_rows == 2){
    if(type(0)==2 && std::abs(df(0,1) - df(1,1)) >= cf && df(1,0) - df(0,0) < s){
      // if only 1 max and 1 min and: min is prior to max, difference in value > cf and they occur within the same year (s)
      // then return the sos and pog
      return List::create(
        Named("df") = df,
        Named("type") = type);
    }else{
      return List::create(
        Named("df") = arma::mat(1,2).fill(NA_REAL),
        Named("type") = NA_REAL);
    }
  }else if(df.n_rows <= 2){ // else return NAs
    return List::create(
      Named("df") = arma::mat(1,2).fill(NA_REAL),
      Named("type") = NA_REAL);
  }else{
    arma::mat df_pog = df.rows(find(type==1));
    arma::mat df_sos = df.rows(find(type==2));
    
    // if sos, pog, and sos identified
    if(df_pog.n_rows == 1){
      if(std::abs(df_sos(0,0) - df_sos(1,0)) < s && std::abs(df_sos(0,1) - df_pog(0,1)) >= cf){
        // if conditions are satisfied, return first sos and pog
        df.shed_row(2);
        type.shed_row(2);
        return List::create(
          Named("df") = df,
          Named("type") = type);
      }else{
        return List::create(
          Named("df") = arma::mat(1,2).fill(NA_REAL),
          Named("type") = NA_REAL);
      }
    }else{
      List t = checkd_min_new(df, type, d); // check pog's for L condition, need to POGs
      arma::mat df2 = t[0]; // remaining times and values
      arma::vec type2 = t[1]; // remaining types
      if(R_IsNA(type2(0)) || type2.n_elem < 2){ // if returned NA or less than 2
        return List::create(
          Named("df") = arma::mat(1,2).fill(NA_REAL),
          Named("type") = NA_REAL);
      }else if(type2.n_elem==2){ // if only two remain
        if(type2(0)==2 && std::abs(df2(0,1) - df2(1,1)) >= cf && df2(1,0) - df2(0,0) < s && std::abs(df2(0,0) - df2(1,0)) > d){
          // if satisfy conditions
          return t;
        }else{
          return List::create(
            Named("df") = arma::mat(1,2).fill(NA_REAL),
            Named("type") = NA_REAL);
        }
      }else if(type2.n_elem==3 && type2(0)==1){ // if ony 3 remaining and first is POG (need two sos for checkd_max_new)
        if(std::abs(df_sos(1,0) - df_sos(2,0)) < s && std::abs(df_sos(1,1) - df_pog(2,1)) >= cf){
          // if conditions are satisfied, return first sos and pog
          df2.shed_row(0);
          type2.shed_row(0);
          return List::create(
            Named("df") = df2,
            Named("type") = type2);
        }else{
          return List::create(
            Named("df") = arma::mat(1,2).fill(NA_REAL),
            Named("type") = NA_REAL);
        }
      }else{
        List t2 = checkd_max_new(df2, type2, d); // check pog's for L condition, need at least two sos
        arma::vec d2 = t2[1];
        if(R_IsNA(d2(0))){ // checkd_max_new return NA if less than 2
          return List::create(
            Named("df") = arma::mat(1,2).fill(NA_REAL),
            Named("type") = NA_REAL);
        }else{
          List t3 = check_2sos(t2[0], t2[1], cf); // checks when prev pog and sos, and sos and next pog are all less than cf -- something like that
          arma::vec d3 = t3[1];
          if(R_IsNA(d3(0))){
            return List::create(
              Named("df") = arma::mat(1,2).fill(NA_REAL),
              Named("type") = NA_REAL);
          }else{
            List t4 = check_2pog(t3[0], t3[1], cf); // checks the same with prev sos and pog, and pog and next sos -- something like that
            arma::vec d4 = t3[1];
            if(R_IsNA(d4(0))){
              return List::create(
                Named("df") = arma::mat(1,2).fill(NA_REAL),
                Named("type") = NA_REAL);
            }else{
              List t5 = check_2up(t4[0], t4[1], y, cf); // checks when previous and next pairs are ok but current are not
              arma::mat df3 = t5[0]; // remaining times and values
              arma::vec type3 = t5[1]; // remaining types
              if(R_IsNA(type3(0)) || type3.n_elem < 2){
                return List::create(
                  Named("df") = arma::mat(1,2).fill(NA_REAL),
                  Named("type") = NA_REAL);
              }else if(type3.n_elem==3 && type3(0)==1){ // if 3 remaining and first is POG (need two SOS for checkd_max_alt)
                if(std::abs(df3(1,0) - df3(2,0)) < s && std::abs(df3(1,1) - df3(2,1)) >= cf){
                  // if conditions are satisfied, return first sos and pog
                  df3.shed_row(0);
                  type3.shed_row(0);
                  return List::create(
                    Named("df") = df3,
                    Named("type") = type3);
                }else{
                  return List::create(
                    Named("df") = arma::mat(1,2).fill(NA_REAL),
                    Named("type") = NA_REAL);
                }
              }else if(type3.n_elem > 2){ // if more than 2 elements
                List t6 = checkd_max_alt(df3, type3, d); // do one last check for L duration
                arma::mat df4 = t6[0]; // remaining times and values
                arma::vec type4 = t6[1]; // remaining types
                if(R_IsNA(type4(0)) || df4.n_rows < 2){
                  return List::create(
                    Named("df") = arma::mat(1,2).fill(NA_REAL),
                    Named("type") = NA_REAL);
                }else if(type4.n_elem==2){
                  if(type4(0)==2 && std::abs(df4(0,1) - df4(1,1)) >= cf && df4(1,0) - df4(0,0) < s && std::abs(df4(0,0) - df4(1,0)) > d){
                    return t6;
                  }else{
                    return List::create(
                      Named("df") = arma::mat(1,2).fill(NA_REAL),
                      Named("type") = NA_REAL);
                  }
                }else{
                  return t6;
                }
              }else{
                if(type3(0)==2 && std::abs(df3(0,1) - df3(1,1)) >= cf && df3(1,0) - df3(0,0) < s && std::abs(df3(0,0) - df3(1,0)) > d){
                  return t5;
                }else{
                  return List::create(
                    Named("df") = arma::mat(1,2).fill(NA_REAL),
                    Named("type") = NA_REAL);
                }
              }
            }
          }
        }
      }
    }
  }
}

// [[Rcpp::export]]
arma::mat add_year_cpp_fixed(List rule_out, double shift, int year_start, int year_end, int s){
  arma::mat df = rule_out[0];
  arma::vec type = rule_out[1];
  
  arma::vec Time = df.col(0);
  arma::vec Year(df.n_rows);
  int n_years = year_end - year_start + 1;
  arma::mat df_f(df.n_rows, 4);
  
  if(R_IsNA(type(0))){
    df_f(span(0, df.n_rows-1), span(0, 1)) = df;
    df_f.col(2) = type;
    df_f.col(3) = NA_REAL;
  }else{
    for(int i=0; i < df.n_rows; i++){
      for(int j=0; j < n_years; j++){
        if(Time(i) <= (s+shift)){
          Year(i) = year_start;
          j = n_years;
        }else if(Time(i) > (s*j + shift) && Time(i) <= (s*(j+1) + shift)){
          Year(i) = year_start + j;
          j = n_years;
        }
      }
    }
    df_f(span(0, df.n_rows-1), span(0, 1)) = df;
    df_f.col(2) = type;
    df_f.col(3) = Year;
  }
  return(df_f);
}

// [[Rcpp::export]]
arma::vec seas_per_year_cpp(arma::vec years, arma::vec type, int year_start, int year_end){
  arma::vec years_pog = years(find(type==1));
  int n_year = year_end - year_start + 1;
  arma::vec n_per_year(n_year);
  arma::vec it;
  for(int j = 0; j < n_year; j++){
    it = years_pog(find(years_pog==(year_start + j)));
    n_per_year(j) = it.n_elem;
  }
  return n_per_year;
}

// [[Rcpp::export]]
double quantile2(arma::vec x, double q, bool upper = false) {
  arma::vec y = sort(x);
  int it = floor(x.n_elem*q);

  if(upper==true){
    return y(x.n_elem-it);
  }else{
    return y(it-1);
  }
}


// [[Rcpp::export]]
List seas_per_year_all_fixed(List pts, double n_iter, int year_start, int year_end, double cut){
   int n_year = year_end - year_start + 1;
   arma::mat all_nyear(n_iter, n_year);
   arma::mat df_i;
   arma::mat probs(n_year,4);
   arma::vec n0(n_year), n1(n_year), n2(n_year), n3(n_year);
   n0.zeros();
   n1.zeros();
   n2.zeros();
   n3.zeros();
   for(int i=0; i < n_iter; i++){
     df_i = as<arma::mat>(pts[i]);
     if(R_IsNA(df_i(0,0))){
       all_nyear.row(i).zeros();
     }else{
       all_nyear.row(i) = seas_per_year_cpp(df_i.col(3), df_i.col(2), year_start, year_end).t();
     }
     for(int j=0; j < n_year; j++){
       if(all_nyear(i,j) == 0){
         n0(j)++;
       }else if(all_nyear(i,j) == 1){
         n1(j)++;
       }else if(all_nyear(i,j) == 2){
         n2(j)++;
       }else if(all_nyear(i,j) == 3){
         n3(j)++;
       }
     }
   }
   probs.col(0) = n0/n_iter;
   probs.col(1) = n1/n_iter;
   probs.col(2) = n2/n_iter;
   probs.col(3) = n3/n_iter;

   arma::vec seas(n_year);
   seas.zeros();

   for(int j=0; j < n_year; j++){
     if(probs(j,3) > cut){
       seas(j) = 3;
     }else if(probs(j,2) > cut){
       seas(j) = 2;
     }else if(probs(j,1) > cut){
       seas(j) = 1;
     }else if(probs(j,0) > cut){
       seas(j) = 0;
     }else{
       seas(j) = armaWhichMax(probs.row(j).t());
     }
   }

   arma::umat ks(n_iter, n_year);
   for(int j = 1; j < n_year; j++){
     if(seas(j)==0){
       ks.col(j) = all_nyear.col(j)==seas(j);
     }else{
       ks.col(j) = all_nyear.col(j)==seas(j) && all_nyear.col(j-1) != 0;
     }
   }

   ks.col(0) = all_nyear.col(0)==seas(0);

   return List::create(
     Named("seas_year") = all_nyear,
     Named("probs") = probs,
     Named("nseas_year") = seas,
     Named("ks") = ks);
 }

// [[Rcpp::export]]
List pts_to_intervals_cpp_fixed(List pts, List seas_output, int year_start, int year_end, double alpha){

  arma::mat all_nyear = seas_output[0];
  arma::mat probs = seas_output[1];
  arma::vec seas = seas_output[2];
  arma::umat ks = seas_output[3];
  int n_year = seas.n_elem;

  List summs(n_year);
  List df1s(2);

  for(int j=1; j < n_year; j++){
    IntegerVector l = as<IntegerVector>(wrap(find(ks.col(j)==1)));
    List kpts_sub = pts[l];
    if(seas(j)==0){
      summs[j] = arma::mat(1,5).zeros();
    }else{
      arma::mat dfl(2*seas(j)*l.length(),5);
      arma::mat df_sum(2*seas(j), 5);
      for(int i=0; i < l.length(); i++){
        arma::mat k2 = kpts_sub[i];
        arma::uvec f = find(k2.col(2)==1 && k2.col(3)==year_start+j); //type == POG for a given year year_start + j
        dfl.submat(2*seas(j)*i,0,2*seas(j)*(i+1)-1,3) = k2.submat(f(0)-1, 0, f(0) + 2*seas(j)-2, 3);
        dfl.submat(2*seas(j)*i,4,2*seas(j)*(i+1)-1,4) = linspace<arma::vec>(1, 2*seas(j), 2*seas(j));
      }
      df_sum.col(0) = dfl.submat(0,2,2*seas(j)-1,2);

      for(int i=1; i<=2*seas(j); i++){
        arma::mat df1 = dfl.rows(find(dfl.col(4)==i));
        df_sum(i-1,1) = quantile2(df1.col(0), alpha/2.0);
        df_sum(i-1,2) = quantile2(df1.col(0), 0.5);
        df_sum(i-1,3) = quantile2(df1.col(0), alpha/2.0, true);
        df_sum(i-1,4) = quantile2(df1.col(1), 0.5);
      }
      summs[j] = df_sum;
    }
  }

  // first year
  IntegerVector l = as<IntegerVector>(wrap(find(ks.col(0)==1)));
  List kpts_sub = pts[l];
  int p = 2*seas(0)-1;
  if(seas(0)==0){
    summs[0] = arma::mat(1,5).zeros();
  }else{
    arma::mat dfl(p*l.length(),5);
    arma::mat df_sum(p, 5);
    for(int i=0; i < l.length(); i++){
      arma::mat k2 = kpts_sub[i];
      arma::uvec f = find(k2.col(2)==1 && k2.col(3)==year_start);
      dfl.submat(p*i,0,p*(i+1)-1,3) = k2.submat(f(0), 0, f(0) + p - 1, 3);
      dfl.submat(p*i,4,p*(i+1)-1,4) = linspace<arma::vec>(1, p, p);
    }
    df_sum.col(0) = dfl.submat(0,2,p-1,2);

    for(int i=1; i<=p; i++){
      arma::mat df1 = dfl.rows(find(dfl.col(4)==i));
      df_sum(i-1,1) = quantile2(df1.col(0), alpha/2.0);
      df_sum(i-1,2) = quantile2(df1.col(0), 0.5);
      df_sum(i-1,3) = quantile2(df1.col(0), alpha/2.0, true);
      df_sum(i-1,4) = quantile2(df1.col(1), 0.5);
    }
    summs[0] = df_sum;
  }

  arma::vec nused(n_year);
  for(int i=0; i < n_year; i++){
    nused(i) = sum(ks.col(i));
  }

  return List::create(
    // Named("seas_year") = all_nyear,
    Named("probs") = probs,
    Named("nseas_year") = seas,
    Named("intervals")=summs,
    Named("nused")=nused);
}

// [[Rcpp::export]]
List rule_to_intervals_cpp_fixed(List pts, double n_iter, int year_start, int year_end, double cut, double alpha){
  List seas_out = seas_per_year_all_fixed(pts, n_iter, year_start, year_end, cut);
  List ints = pts_to_intervals_cpp_fixed(pts, seas_out, year_start, year_end, alpha);
  return ints;
}

// [[Rcpp::export]]
List mcmc_to_intervals_cpp_fixed(arma::mat St, double cf, int d, double shift, int s, int year_start, int year_end, double cut, double alpha, bool return_points = false){
  List pts(St.n_cols);
  for(int i=0; i < St.n_cols; i++){
    List temp = rule_cpp_no_errors(St.col(i), cf, d, s);
    pts[i] = add_year_cpp_fixed(temp, shift, year_start, year_end, s);
  }
  List seas_out = seas_per_year_all_fixed(pts, St.n_cols, year_start, year_end, cut);
  List ints = pts_to_intervals_cpp_fixed(pts, seas_out, year_start, year_end, alpha);
  if(return_points==true){
    return List::create(
      Named("ints") = ints,
      Named("pts") = pts);
  }else{
    return List::create(
      Named("ints") = ints);
  }
}


