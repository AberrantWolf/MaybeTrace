use std::ops;

#[derive(Copy, Clone)]
struct Vec3D {x:f32, y:f32, z:f32}
impl Vec3D {
    fn new(ix: f32, iy: f32, iz: f32) -> Vec3D {
        Vec3D { x: ix, y: iy, z: iz }
    }

    fn dot(self, other: Vec3D) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn cross(&self, other: Vec3D) -> Vec3D {
        Vec3D {
            x: self.y * other.z - self.z * other.y,
            y: -self.x * other.z - self.z * other.x,
            z: self.x * other.y - self.y * other.x
        }
    }

    fn squared_length(&self) -> f32 {
        self.x*self.x + self.y*self.y + self.z*self.z
    }

    fn length(&self) -> f32 {
        self.squared_length().sqrt()
    }

    fn as_unit(&self) -> Vec3D {
        let mag = self.length();
        Vec3D {
            x: self.x / mag,
            y: self.y / mag,
            z: self.z / mag
        }
    }

    fn normalize(&mut self) {
        let len = self.length();
        self.x /= len;
        self.y /= len;
        self.z /= len;
    }
}

impl ops::Add<Vec3D> for Vec3D {
    type Output = Vec3D;

    fn add(self, _rhs: Vec3D) -> Vec3D {
        Vec3D {
            x: self.x + _rhs.x,
            y: self.y + _rhs.y,
            z: self.z + _rhs.z
        }
    }
}

impl ops::Sub<Vec3D> for Vec3D {
    type Output = Vec3D;

    fn sub(self, _rhs: Vec3D) -> Vec3D {
        Vec3D {
            x: self.x - _rhs.x,
            y: self.y - _rhs.y,
            z: self.z - _rhs.z,
        }
    }
}

impl ops::Mul<Vec3D> for Vec3D {
    type Output = Vec3D;

    fn mul(self, _rhs: Vec3D) -> Vec3D {
        Vec3D {
            x: self.x * _rhs.x,
            y: self.y * _rhs.y,
            z: self.z * _rhs.z,
        }
    }
}

impl ops::Div<Vec3D> for Vec3D {
    type Output = Vec3D;

    fn div(self, _rhs: Vec3D) -> Vec3D {
        Vec3D {
            x: self.x / _rhs.x,
            y: self.y / _rhs.y,
            z: self.z / _rhs.z,
        }
    }
}

impl ops::Mul<f32> for Vec3D {
    type Output = Vec3D;

    fn mul(self, _rhs: f32) -> Vec3D {
        Vec3D {
            x: self.x * _rhs,
            y: self.y * _rhs,
            z: self.z * _rhs,
        }
    }
}

impl ops::Div<f32> for Vec3D {
    type Output = Vec3D;

    fn div(self, _rhs: f32) -> Vec3D {
        Vec3D {
            x: self.x / _rhs,
            y: self.y / _rhs,
            z: self.z / _rhs,
        }
    }
}

impl ops::Mul<Vec3D> for f32 {
    type Output = Vec3D;

    fn mul(self, _rhs: Vec3D) -> Vec3D {
        _rhs * self
    }
}

impl ops::Div<Vec3D> for f32 {
    type Output = Vec3D;

    fn div(self, _rhs: Vec3D) -> Vec3D {
        _rhs * self
    }
}

struct Ray {
    origin: Vec3D,
    direction: Vec3D,
}
impl Ray {
    fn new(orig: Vec3D, dir: Vec3D) -> Ray {
        Ray {origin:orig, direction:dir}
    }

    fn point_at_dist(self, dist:f32) -> Vec3D {
        self.origin + self.direction * dist
    }
}

struct Dimension {x:usize, y:usize}

struct Image {
    dims: Dimension,
    data: Vec<Vec3D>,
}

impl Image {
    pub fn new(dim: Dimension) -> Image {
        let mut result = Image {
            dims: dim,
            data: Vec::new()
        };

        for _y in 0..result.dims.y {
            for _x in 0..result.dims.x {
                result.data.push(Vec3D{
                    x: 0f32,
                    y: 0f32,
                    z: 0f32,
                })
            }
        }

        return result;
    }
}

struct HitRecord {
    t:f32,
    p:Vec3D,
    normal:Vec3D,
}
trait Hitable {
    fn hit(ray:&Ray, t_min:f32, t_max:f32, hit_record:&mut HitRecord) -> bool;
}

// TODO: create Sphere class and implement Hitable on it
struct Sphere {
    center:Vec3D,
    radius:f32,
}

fn hit_sphere(center:Vec3D, radius:f32, ray:&Ray) -> f32 {
    let oc = ray.origin - center;
    let a = ray.direction.dot(ray.direction);
    let b = 2.0 * oc.dot(ray.direction);
    let c = oc.dot(oc) - radius * radius;
    let descriminant = b*b - 4.0*a*c;
    if descriminant < 0.0 {
        return -1.0;
    }

    (-b - descriminant.sqrt()) / (2.0 * a)
}

fn color(ray:Ray) -> Vec3D {
    let t = hit_sphere(Vec3D::new(0.0, 0.0, -1.0), 0.5, &ray);

    if t > 0.0 {
        let n = (ray.point_at_dist(t) - Vec3D::new(0.0, 0.0, -1.0)).as_unit();
        return 0.5 * Vec3D::new(n.x + 1.0, n.y + 1.0, n.z + 1.0);
    }

    let unit_dir:Vec3D = ray.direction.as_unit();
    let t:f32 = 0.5 * (unit_dir.y + 1.0);
    (1.0 - t) * Vec3D::new(1.0, 1.0, 1.0) + t * Vec3D::new(0.5, 0.7, 1.0)
}

fn main() {
    let mut img = Image::new(Dimension{x:200, y:100});

    let lower_left_corner = Vec3D::new(-2.0, -1.0, -1.0);
    let horizontal = Vec3D::new(4.0, 0.0, 0.0);
    let vertical = Vec3D::new(0.0, 2.0, 0.0);
    let origin = Vec3D::new(0.0, 0.0, 0.0);

    println!("P3");
    println!("{0} {1}", img.dims.x, img.dims.y);
    println!("255");
    for y in 0..img.dims.y {
        for x in 0..img.dims.x {
            let u = x as f32 / img.dims.x as f32;
            let v = 1f32 - y as f32 / img.dims.y as f32;

            let r = Ray::new(origin, lower_left_corner + u*horizontal + v*vertical);
            let c = color(r);

            println!("{0} {1} {2}",
                     (c.x * 255.99) as i32,
                     (c.y * 255.99) as i32,
                     (c.z * 255.99) as i32
            );
            &img.data.push(c);
        }
    }
}
