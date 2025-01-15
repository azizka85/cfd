use std::fs;
use std::io::{self, Write};

const U: f64 = 1.;

const H: f64 = 1.;

const DIR: &str = "data/1d/couette_flow/steady_state";

fn main() {        
    if H < 0. {
        println!("H should be >= 0, but it is {H}");

        return;
    }

    let dy = 0.01;

    if dy <= 0. {
        println!("dy should be > 0, but it is {dy}");

        return;
    }

    let ny = (H/dy).ceil() as usize + 1;

    let res = set_environment(dy);

    println!("Set Environment: {:?}", res);

    let path = match res {
        Err(_) => return,
        Ok(path) => path
    };
    
    let mut u = vec![0.; ny];

    set_boundary_condition(&mut u);

    let [l, c, r] = construct_slae(ny);

    u = match slae::direct::tridiagonal::solve(&l, &c, &r, &mut u) {
        Err(err) => {
            println!("Error in solving tridiagonal matrix: {err}");

            return;
        },
        Ok(u) => u
    };


    if let Err(err) = write_file(&u, dy, &path) {
        println!("Error in writing file: {:?}", err);

        return;
    }
}

fn set_environment(dy: f64) -> io::Result<String> {
    let path = format!("{DIR}/U={U}, H={H}/dy={dy}");

    fs::create_dir_all(path.clone())?;

    Ok(path)
}

fn write_file(u: &Vec<f64>, dy: f64, path: &str) -> io::Result<()> {
    let mut file = fs::File::create(format!("{path}/data.vtk"))?;

    writeln!(file, "# vtk DataFile Version 3.0")?;
    writeln!(file, "TIME 0")?;
    writeln!(file, "ASCII")?;
    writeln!(file, "DATASET STRUCTURED_GRID")?;

    let ny = u.len();

    writeln!(file, "DIMENSIONS 1 {ny} 1")?;
    writeln!(file, "POINTS {ny} float")?;

    for i in 0..ny {
        let y = dy * i as f64;

        writeln!(file, "0.0 {y:.3} 0.0")?;
    }    

    writeln!(file, "POINT_DATA {}", ny)?;
    writeln!(file, "SCALARS u float")?;
    writeln!(file, "LOOKUP_TABLE default")?;

    for i in 0..ny {
        writeln!(file, "{:.3}", u[i])?;
    } 

    Ok(())
}

fn set_boundary_condition(u: &mut Vec<f64>) {
    let n = u.len();

    if n > 0 {
        u[0] = 0.;
        u[n-1] = U;
    }
}

fn construct_slae(ny: usize) -> [Vec<f64>; 3] {
    let mut l = vec![0.; ny];
    let mut c = vec![0.; ny];
    let mut r = vec![0.; ny];

    if ny > 0 {
        c[0] = 1.;
        c[ny-1] = 1.;

        for i in 1..ny-1 {
            l[i] = 1.;
            c[i] = -2.;
            r[i] = 1.;
        }
    }
    
    [l, c, r]
}
