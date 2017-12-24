use std::ops;
use std::f32;

#[derive(Copy, Clone)]
struct Vec3D {x:f32, y:f32, z:f32}
impl Vec3D {
    fn new(ix: f32, iy: f32, iz: f32) -> Vec3D {
        Vec3D { x: ix, y: iy, z: iz }
    }

    fn dot(&self, other: Vec3D) -> f32 {
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

    fn point_at_dist(&self, dist:f32) -> Vec3D {
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
impl HitRecord {
    fn new() -> HitRecord {
        HitRecord {
            t: 0.0,
            p: Vec3D { x: 0.0, y: 0.0, z: 0.0 },
            normal: Vec3D { x: 0.0, y: 0.0, z: 0.0 },
        }
    }
}

trait Hitable {
    fn hit(&self, ray:&Ray, t_min:f32, t_max:f32, hit_record:&mut HitRecord) -> bool;
}

struct Sphere {
    center:Vec3D,
    radius:f32,
}
impl Sphere {
    fn new(ctr:Vec3D, rad:f32) -> Sphere {
        Sphere {center:ctr, radius:rad}
    }
}

impl Hitable for Sphere {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32, hit_record: &mut HitRecord) -> bool {
        let oc = ray.origin - self.center;
        let a = ray.direction.dot(ray.direction);
        let b = oc.dot(ray.direction);
        let c = oc.dot(oc) - self.radius * self.radius;

        let descriminant = b*b - a*c;
        if descriminant > 0.0 {
            let t_neg = (-b - descriminant.sqrt()) / a;
            let t_pos = (-b + descriminant.sqrt()) / a;

            if t_neg < t_max && t_neg > t_min {
                hit_record.t = t_neg;
                hit_record.p = ray.point_at_dist(hit_record.t);
                hit_record.normal = (hit_record.p - self.center) / self.radius;
                return true
            }
            else if t_pos < t_max && t_pos > t_min {
                hit_record.t = t_pos;
                hit_record.p = ray.point_at_dist(hit_record.t);
                hit_record.normal = (hit_record.p - self.center) / self.radius;
                return true
            }
        }
        return false
    }
}

struct HitableList {
    list:Vec<Box<Hitable>>,
}

impl HitableList {
    fn new() -> HitableList {
        HitableList {
            list: Vec::new(),
        }
    }
}

impl Hitable for HitableList {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32, hit_record: &mut HitRecord) -> bool {
        let mut hit_anything = false;
        let mut closest_so_far = t_max;

        for h in &self.list {
            let mut temp_record = HitRecord::new();
            if h.hit(ray, t_min, closest_so_far, &mut temp_record) {
                hit_anything = true;
                closest_so_far = temp_record.t;
                *hit_record = temp_record;
            }
        }

        return hit_anything;
    }
}

fn color(ray:Ray, world:&Hitable) -> Vec3D {
    let mut rec = HitRecord::new();
    if world.hit(&ray, 0.0, f32::MAX, &mut rec) {
        return 0.5 * Vec3D::new(rec.normal.x + 1.0, rec.normal.y + 1.0, rec.normal.z + 1.0);
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

    let mut world = HitableList::new();
    world.list.push(Box::new(Sphere::new(Vec3D::new(0.0, 0.0, -1.0), 0.5)));
    world.list.push(Box::new(Sphere::new(Vec3D::new(0.0, -100.5, -1.0), 100.0)));

    println!("P3");
    println!("{0} {1}", img.dims.x, img.dims.y);
    println!("255");
    for y in 0..img.dims.y {
        for x in 0..img.dims.x {
            let u = x as f32 / img.dims.x as f32;
            let v = 1f32 - y as f32 / img.dims.y as f32;

            let r = Ray::new(origin, lower_left_corner + u*horizontal + v*vertical);
            let c = color(r, &world);

            println!("{0} {1} {2}",
                     (c.x * 255.99) as i32,
                     (c.y * 255.99) as i32,
                     (c.z * 255.99) as i32
            );
            &img.data.push(c);
        }
    }
}
