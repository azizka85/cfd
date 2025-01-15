pub fn solve(l: &Vec<f64>, c: &Vec<f64>, r: &Vec<f64>, d: &mut Vec<f64>) -> Result<Vec<f64>, String> {
    let n = c.len();    

    if l.len() != n {
        return Err(
            format!("The lengths of the l and c should be equal, but now the length of the l is {} and c is {n}", l.len())
        );
    }

    if r.len() != n {
        return Err(
            format!("The lengths of the r and c should be equal, but now the length of the r is {} and c is {n}", r.len())
        );
    }

    if d.len() != n {
        return Err(
            format!("The lengths of the d and c should be equal, but now the length of the d is {} and c is {n}", d.len())
        );
    }

    let mut u = vec![0.; n];

    if n > 0 {
        u[0] = r[0] / c[0];
        d[0] = d[0] / c[0];                

        for i in 1..n {
            let c1 = c[i] - l[i]*u[i-1];
    
            u[i] = r[i] / c1;
            d[i] = (d[i] - l[i]*d[i-1]) / c1;
        }

        u[n-1] = d[n-1];

        for i in (0..n-1).rev() {
            u[i] = d[i] - u[i]*u[i+1];
        }
    }    

    Ok(u)
}
