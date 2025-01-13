use std::fs;
use std::io::{self, Write};

const NU: f64 = 1.;
const U: f64 = 1.;

const H: f64 = 1.;
const T: f64 = 0.3;

const OTS: f64 = 0.003;

const DIR: &str = "data/1d/couette_flow/explicit_euler";

fn main() {        
    let dy = 0.01;
    let ny = (H/dy).ceil() as usize + 1;

    let r = 0.3; 
    let dt = r*dy*dy/NU;   

    let res = set_environment(dy, r);

    println!("Set Environment: {:?}", res);

    let path = match res {
        Err(_) => return,
        Ok(path) => path
    };
    
    let mut u = vec![0.; ny];
    let mut u1 = vec![0.; ny];

    set_initial_condition(&mut u);
    set_boundary_condition(&mut u);

    let mut t = 0.;
    let mut tn = OTS;

    let mut n = 0;    

    if let Err(err) = write_file(&u, t, dy, n, &path) {
        println!("Error in writing file: n={n}, t={t}: {:?}", err);

        return;
    } else {
        n += 1;
    }
    
    while t <= T {
        set_boundary_condition(&mut u1);

        for i in 1..ny-1 {
            u1[i] = u[i] + r*(u[i+1] - 2.*u[i] + u[i-1]);
        }

        update_data(&mut u, &mut u1);

        t += dt;                        

        if t >= tn {
            if let Err(err) = write_file(&u, t, dy, n, &path) {
                println!("Error in writing file: n={n}, t={t}: {:?}", err);
        
                return;
            } else {
                n += 1;
            }

            tn += OTS;
        }
    } 
}

fn set_environment(dy: f64, r: f64) -> io::Result<String> {
    let path = format!("{DIR}/nu={NU}, U={U}, H={H}/dy={dy}, r={r}");

    fs::create_dir_all(path.clone())?;

    Ok(path)
}

fn write_file(u: &Vec<f64>, t: f64, dy: f64, n: usize, path: &str) -> io::Result<()> {
    let mut file = fs::File::create(format!("{path}/data.{n:03}.vtk"))?;

    writeln!(file, "# vtk DataFile Version 3.0")?;
    writeln!(file, "TIME {t:.3}")?;
    writeln!(file, "ASCII")?;
    writeln!(file, "DATASET STRUCTURED_GRID")?;

    let ny = u.len();

    writeln!(file, "DIMENSIONS 1 {ny} 1")?;
    writeln!(file, "POINTS {ny} float")?;

    for i in 0..ny {
        let y = dy * i as f64;

        writeln!(file, "0.0 {y:.3} 0.0")?;
    }

    writeln!(file, "FIELD FieldData 1")?;
    writeln!(file, "Time 1 1 float")?;
    writeln!(file, "{t:.3}")?;

    writeln!(file, "POINT_DATA {}", ny)?;
    writeln!(file, "SCALARS u float")?;
    writeln!(file, "LOOKUP_TABLE default")?;

    for i in 0..ny {
        writeln!(file, "{:.3}", u[i])?;
    } 

    Ok(())
}

fn set_initial_condition(u: &mut Vec<f64>) {
    let n = u.len();

    for i in 0..n {
        u[i] = 0.;
    }
}

fn set_boundary_condition(u: &mut Vec<f64>) {
    let n = u.len();

    if n > 0 {
        u[0] = 0.;
        u[n-1] = U;
    }
}

fn update_data(u: &mut Vec<f64>, u1: &mut Vec<f64>) {
    let n = u.len();

    for i in 0..n {
        u[i] = u1[i];
        u1[i] = 0.;
    }
}
