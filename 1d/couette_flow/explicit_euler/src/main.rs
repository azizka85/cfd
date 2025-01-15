use std::fs;
use std::io::{self, Write};

const NU: f64 = 1.;
const U: f64 = 1.;

const H: f64 = 1.;
const T: f64 = 0.3;

const OTS: f64 = 0.003;

const DIR: &str = "data/1d/couette_flow/explicit_euler";

fn main() {        
    if NU < 0. {
        println!("NU should be >= 0, but it is {NU}");

        return;
    }

    if H < 0. {
        println!("H should be >= 0, but it is {H}");

        return;
    }

    if OTS <= 0. {
        println!("OTS should be > 0, but it is {OTS}");

        return;
    }

    let dy = 0.01;

    if dy <= 0. {
        println!("dy should be > 0, but it is {dy}");

        return;
    }

    let ny = (H/dy).ceil() as usize + 1;

    let r = 0.3; 

    if r <= 0. {
        println!("r should be > 0, but it is {r}");

        return;
    }

    let dt = r*dy*dy/NU;   

    let res = set_environment(dy, r);

    println!("Set Environment: {:?}", res);

    let path = match res {
        Err(_) => return,
        Ok(path) => path
    };
    
    let mut u = vec![0.; ny];
    let mut u1 = vec![0.; ny];    
    let mut up = vec![0.; ny];

    set_initial_condition(&mut u);
    set_boundary_condition(&mut u);     

    update_data(&mut up, &u);   

    let mut t = 0.;
    let mut tn = OTS;

    let mut n = 1;
    let mut m = 0; 

    let mut statistics = Vec::new();   

    if let Err(err) = write_data(&u, t, dy, m, &path) {
        println!("Error in writing file: m={m}, t={t}: {:?}", err);

        return;
    } else {
        m += 1;
    }
    
    while t <= T {
        set_boundary_condition(&mut u1);

        for i in 1..ny-1 {
            u1[i] = u[i] + r*(u[i+1] - 2.*u[i] + u[i-1]);
        }        

        update_data(&mut u, &u1);
        set_value(&mut u1, 0.);

        t += dt;                        

        if t >= tn {
            if let Err(err) = write_data(&u, t, dy, m, &path) {
                println!("Error in writing file: m={m}, t={t}: {:?}", err);
        
                return;
            } else {
                let max_diff = max_abs_difference(&u, &up);

                println!("Write data in file t={t:.3}, convergence of u={max_diff:.5}");                

                statistics.push((n, tn, max_diff));

                update_data(&mut up, &u);

                m += 1;
            }

            tn += OTS;
        }

        n += 1;
    } 

    let res = write_statistics(&statistics, &path);

    println!("Write Statistics: {:?}", res);
}

fn set_environment(dy: f64, r: f64) -> io::Result<String> {
    let path = format!("{DIR}/nu={NU}, U={U}, H={H}/dy={dy}, r={r}");

    fs::create_dir_all(path.clone())?;

    Ok(path)
}

fn write_data(u: &Vec<f64>, t: f64, dy: f64, m: usize, path: &str) -> io::Result<()> {
    let mut file = fs::File::create(format!("{path}/data.{m:03}.vtk"))?;

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

fn write_statistics(statistics: &Vec<(usize, f64, f64)>, path: &str) -> io::Result<()> {
    let path = format!("{path}/statistics");
    
    fs::create_dir_all(path.clone())?;

    let mut writer = csv::Writer::from_path(
        format!("{path}/convergence.csv")
    )?;
    
    writer.write_record(["n", "t", "max_diff"])?;

    for s in statistics {
        writer.write_record([
            format!("{}", s.0),
            format!("{:.3}", s.1),
            format!("{:.5}", s.2)
        ])?;
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

fn update_data(u: &mut Vec<f64>, u1: &Vec<f64>) {
    let n = u.len();

    for i in 0..n {
        u[i] = u1[i];
    }
}

fn set_value(u: &mut Vec<f64>, value: f64) {
    let n = u.len();

    for i in 0..n {
        u[i] = value;
    }
}

fn max_abs_difference(u: &Vec<f64>, u1: &Vec<f64>) -> f64 {
    let n = u.len();

    let mut max_diff: f64 = 0.;

    for i in 0..n {
        max_diff = max_diff.max((u[i] - u1[i]).abs())
    }

    max_diff
}
