use std::fs::read_to_string;

pub fn input_mesh(mesh_file: &str) {
    // Read all lines
    let lines_all =
        match read_to_string(mesh_file) {
            Ok(data) => data,
            Err(err) => {panic!("Could not open mesh file: {}", err);}
        };
    let lines: Vec<&str> = lines_all.trim().split("\n").collect();

    // Read header
    let mut items: Vec<&str> = lines[0].split_whitespace().collect();
    let nnode: usize        = items[0].parse().unwrap();
    let nelem: usize        = items[1].parse().unwrap();
    let nmaterial: usize    = items[2].parse().unwrap();
    let dof: usize          = items[3].parse().unwrap();
    // println!("{} {} {} {}",nnode,nelem,nmaterial,dof);

    // Read node
    let mut irec = 1;
    for inode in 0..nnode {
        items = lines[irec+inode].split_whitespace().collect();

        let id: usize = items[0].parse().unwrap();
        let mut xyz: Vec<f64> = Vec::new();
        for i in 1..3 {
            xyz.push(items[i].parse().unwrap());
        }
        let mut freedom: Vec<usize> = Vec::new();
        for i in 3..items.len() {
            freedom.push(items[i].parse().unwrap());
        }

        let node = node::Node {id: &id};

        println!("{} {:?} {:?}",id,xyz,freedom);
    }

    // Read element
    irec += nnode;
    for ielem in 0..nelem {
        items = lines[irec+ielem].split_whitespace().collect();

        let id: usize = items[0].parse().unwrap();
        let style: &str = items[1];
        let material_id: isize = items[2].parse().unwrap();
        let mut inode: Vec<usize> = Vec::new();
        for i in 3..items.len() {
            inode.push(items[i].parse().unwrap());
        }

        println!("{} {} {} {:?}",id,style,material_id,inode);
    }

    // Read material
    irec += nelem;
    for imaterial in 0..nmaterial {
        items = lines[irec+imaterial].split_whitespace().collect();

        let id: usize = items[0].parse().unwrap();
        let style: &str = items[1];
        let mut param: Vec<f64> = Vec::new();
        for i in 2..items.len() {
            param.push(items[i].parse().unwrap());
        }

        println!("{} {} {:?}",id,style,param);
    }

    return;

}
