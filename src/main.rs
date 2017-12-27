use std::ops;
use std::f32;
use std::rc::Rc;

extern crate rand;
use rand::random;

#[derive(Copy, Clone)]
struct Vec3D {x:f32, y:f32, z:f32}
impl Vec3D {
    fn new(ix: f32, iy: f32, iz: f32) -> Vec3D {
        Vec3D { x: ix, y: iy, z: iz }
    }

    fn dot(&self, other: &Vec3D) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn cross(&self, other: &Vec3D) -> Vec3D {
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

impl std::ops::AddAssign for Vec3D {
    fn add_assign(&mut self, other:Vec3D) {
        *self = Vec3D {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        };
    }
}

impl std::ops::SubAssign for Vec3D {
    fn sub_assign(&mut self, other:Vec3D) {
        *self = Vec3D {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        };
    }
}

impl std::ops::MulAssign<f32> for Vec3D {
    fn mul_assign(&mut self, other:f32) {
        *self = Vec3D {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        };
    }
}

impl std::ops::DivAssign<f32> for Vec3D {
    fn div_assign(&mut self, other:f32) {
        *self = Vec3D {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        };
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

#[derive(Copy, Clone)]
struct Ray {
    origin: Vec3D,
    direction: Vec3D,
}
impl Ray {
    fn new(orig: Vec3D, dir: Vec3D) -> Ray {
        Ray {origin:orig, direction:dir}
    }
    fn forward() -> Ray {
        Ray {
            origin: Vec3D { x: 0.0, y: 0.0, z: 0.0 },
            direction: Vec3D { x: 0.0, y: 0.0, z: -1.0 },
        }
    }

    fn point_at_dist(&self, dist:f32) -> Vec3D {
        self.origin + self.direction * dist
    }
}

#[derive(Copy, Clone)]
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

        result
    }
}

struct HitRecord {
    t:f32,
    p:Vec3D,
    normal:Vec3D,
    material:Option<Rc<Material>>,
}
impl HitRecord {
    fn new() -> HitRecord {
        HitRecord {
            t: 0.0,
            p: Vec3D { x: 0.0, y: 0.0, z: 0.0 },
            normal: Vec3D { x: 0.0, y: 0.0, z: 0.0 },
            material: None,
        }
    }
}

trait Hitable {
    fn hit(&self, ray:&Ray, t_min:f32, t_max:f32) -> Option<HitRecord>;
}

struct Sphere<> {
    center:Vec3D,
    radius:f32,
    material:Rc<Material>,
}
impl Sphere {
    fn new(ctr:Vec3D, rad:f32, mat:Rc<Material>) -> Sphere {
        Sphere {
            center:ctr,
            radius:rad,
            material:mat,
        }
    }
}

impl Hitable for Sphere {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let oc = ray.origin - self.center;
        let a = ray.direction.dot(&ray.direction);
        let b = oc.dot(&ray.direction);
        let c = oc.dot(&oc) - self.radius * self.radius;

        let descriminant = b*b - a*c;
        if descriminant > 0.0 {
            let t_neg = (-b - descriminant.sqrt()) / a;
            let t_pos = (-b + descriminant.sqrt()) / a;

            if t_neg < t_max && t_neg > t_min {
                let rp = ray.point_at_dist(t_neg);
                return Some(
                    HitRecord {
                        t: t_neg,
                        p: rp,
                        normal: (rp - self.center) / self.radius,
                        material: Some(Rc::clone(&self.material)),
                    })
            } else if t_pos < t_max && t_pos > t_min {
                let rp = ray.point_at_dist(t_pos);
                return Some(
                    HitRecord {
                        t: t_pos,
                        p: rp,
                        normal: (rp - self.center) / self.radius,
                        material: Some(Rc::clone(&self.material)),
                    })
            }
        }

        None
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
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let mut hit_anything: Option<HitRecord> = None;
        let mut closest_so_far = t_max;

        for h in &self.list {
            if let Some(hr) = h.hit(ray, t_min, closest_so_far) {
                closest_so_far = hr.t;
                drop(hit_anything);
                hit_anything = Some(hr);
            }
        }

        hit_anything
    }
}

struct Camera {
    origin:Vec3D,
    lower_left_corner:Vec3D,
    horizontal:Vec3D,
    vertical:Vec3D,
}
impl Camera {
    fn new() -> Camera {
        Camera {
            origin:Vec3D{x:0.0, y:0.0, z:0.0},
            lower_left_corner:Vec3D{x:-2.0, y:-1.0, z:-1.0},
            horizontal:Vec3D{x:4.0, y:0.0, z:0.0},
            vertical:Vec3D{x:0.0, y:2.0, z:0.0},
        }
    }

    fn get_ray(&self, u:f32, v:f32) -> Ray {
        Ray {origin:self.origin, direction:self.lower_left_corner + u*self.horizontal + v*self.vertical - self.origin}
    }
}

fn random_in_unit_sphere() -> Vec3D {
    let mut p:Vec3D;
    loop {
        p = 2.0 * Vec3D::new(random::<f32>(), random::<f32>(), random::<f32>()) - Vec3D::new(1.0, 1.0, 1.0);

        if p.squared_length() < 1.0 { break; }
    }

    p
}

fn reflect(v:&Vec3D, n:&Vec3D) -> Vec3D {
    *v - 2.0 * v.dot(n) * *n
}

trait Material {
    fn scatter(&self, r_in:&Ray, rec:&HitRecord, attenuation:&mut Vec3D, scattered:&mut Ray) -> bool;
}

struct LambertianMat {
    albedo: Vec3D,
}
impl Material for LambertianMat {
    fn scatter(&self, _r_in:&Ray, rec:&HitRecord, attenuation:&mut Vec3D, scattered:&mut Ray) -> bool {
        let target = rec.p + rec.normal + random_in_unit_sphere();
        *scattered = Ray { origin: rec.p, direction: target - rec.p };
        *attenuation = self.albedo;

        true
    }
}

struct MetalMat {
    albedo:Vec3D
}
impl Material for MetalMat {
    fn scatter(&self, r_in:&Ray, rec:&HitRecord, attenuation:&mut Vec3D, scattered:&mut Ray) -> bool {
        let reflected = reflect(&r_in.direction.as_unit(), &rec.normal);
        *scattered = Ray { origin: rec.p, direction: reflected };
        *attenuation = self.albedo;

        scattered.direction.dot(&rec.normal) > 0.0
    }
}

fn color(ray:Ray, world:&Hitable, depth:i32) -> Box<Vec3D> {
    if let Some(mut rec) = world.hit(&ray, 0.001, f32::MAX) {
        let mut scattered = Ray::forward();
        let mut attenuation = Vec3D { x: 0.0, y: 0.0, z: 0.0 };

        let mut m = rec.material.clone();
        if let Some(ref mat) = m {
            if depth < 5 && mat.scatter(&ray, &mut rec, &mut attenuation, &mut scattered) {
                Box::new(attenuation * *color(scattered, world, depth + 1))
            } else {
                Box::new(Vec3D { x: 0.0, y: 0.0, z: 0.0 })
            }
        } else {
            Box::new(Vec3D { x: 0.0, y: 0.0, z: 0.0 })
        }
    } else {
        let unit_dir: Vec3D = ray.direction.as_unit();
        let t: f32 = 0.5 * (unit_dir.y + 1.0);
        Box::new((1.0 - t) * Vec3D::new(1.0, 1.0, 1.0) + t * Vec3D::new(0.5, 0.7, 1.0))
    }
}

fn main() {
    let mut img = Image::new(Dimension { x: 200, y: 100 });
    let ns: i32 = 100;

    let cam = Camera::new();

    let mut world = HitableList::new();
    world.list.push(Box::new(Sphere::new(Vec3D::new(0.0, 0.0, -1.0), 0.5, Rc::new(LambertianMat { albedo: Vec3D { x: 0.8, y: 0.3, z: 0.3 } }))));
    world.list.push(Box::new(Sphere::new(Vec3D::new(0.0, -100.5, -1.0), 100.0, Rc::new(LambertianMat { albedo: Vec3D { x: 0.8, y: 0.8, z: 0.0 } }))));

    // metal spheres
    world.list.push(Box::new(Sphere::new(Vec3D::new( 1.0, 0.0, -1.0), 0.5, Rc::new(MetalMat { albedo: Vec3D { x: 0.8, y: 0.6, z: 0.2 } }))));
    world.list.push(Box::new(Sphere::new(Vec3D::new(-1.0, 0.0, -1.0), 0.5, Rc::new(MetalMat { albedo: Vec3D { x: 0.8, y: 0.8, z: 0.8 } }))));

    println!("P3");
    println!("{0} {1}", img.dims.x, img.dims.y);
    println!("255");

    for y in 0..img.dims.y {
        for x in 0..img.dims.x {
            let mut col = Vec3D { x: 0.0, y: 0.0, z: 0.0 };
            for _s in 0..ns {
                let u = (x as f32 + random::<f32>()) / img.dims.x as f32;
                let v = 1f32 - (y as f32 + random::<f32>()) / img.dims.y as f32;
                let r = cam.get_ray(u, v);

                col += *color(r, &world, 0);
            }

            col /= ns as f32;
            // gamma correction
            let gammaa_col = Vec3D { x: col.x.sqrt(), y: col.y.sqrt(), z: col.z.sqrt() };
            println!("{0} {1} {2}",
                     (gammaa_col.x * 255.99) as i32,
                     (gammaa_col.y * 255.99) as i32,
                     (gammaa_col.z * 255.99) as i32
            );
            &img.data.push(gammaa_col);
        }
    }
}
