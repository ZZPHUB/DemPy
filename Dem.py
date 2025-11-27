import math
import numpy as np
import vtk
import os

class Particle:
    """单个粒子类"""
    def __init__(self, pos, vel, radius, mass):
        self.position = np.array(pos, dtype=float)
        self.velocity = np.array(vel, dtype=float)
        self.force = np.zeros(3)
        self.radius = radius
        self.mass = mass
        self.fixed = False  # 是否固定位置

class Cube:
    """正方体类，由多个粒子组成"""
    def __init__(self, center, size, particle_radius, density, velocity=(0, 0, 0)):
        self.center = np.array(center, dtype=float)
        self.size = size
        self.particles = []
        self.velocity = np.array(velocity, dtype=float)
        
        # 计算每个粒子的质量（铝的密度约2700 kg/m³）
        particle_volume = (4/3) * math.pi * particle_radius**3
        particle_mass = density * particle_volume
        
        # 在正方体内创建粒子网格
        particles_per_side = max(3, int(size / (2 * particle_radius)))
        spacing = size / particles_per_side
        
        for i in range(particles_per_side):
            for j in range(particles_per_side):
                for k in range(particles_per_side):
                    # 计算粒子在正方体内的相对位置
                    x = self.center[0] - size/2 + spacing/2 + i * spacing
                    y = self.center[1] - size/2 + spacing/2 + j * spacing
                    z = self.center[2] - size/2 + spacing/2 + k * spacing
                    
                    particle = Particle(
                        [x, y, z],
                        velocity,
                        particle_radius,
                        particle_mass
                    )
                    self.particles.append(particle)

class DEMSimulator:
    """离散元模拟器"""
    def __init__(self):
        self.particles = []
        self.cubes = []
        self.dt = 1e-5  # 时间步长
        self.gravity = np.array([0, 0, -9.81])  # 重力加速度
        
        # 材料参数（铝）
        self.youngs_modulus = 70e9  # 杨氏模量 (Pa)
        self.poissons_ratio = 0.33  # 泊松比
        self.restitution_coef = 0.7  # 恢复系数
        self.friction_coef = 0.3    # 摩擦系数
        
    def add_cube(self, center, size, particle_radius, density=2700, velocity=(0, 0, 0)):
        """添加一个正方体"""
        cube = Cube(center, size, particle_radius, density, velocity)
        self.cubes.append(cube)
        self.particles.extend(cube.particles)
        return cube
    
    def calculate_contact_force(self, p1, p2):
        """计算两个粒子之间的接触力"""
        # 计算相对位置和距离
        r_vec = p2.position - p1.position
        distance = np.linalg.norm(r_vec)
        
        # 如果没有接触，返回零力
        if distance > p1.radius + p2.radius:
            return np.zeros(3)
        
        # 计算法向重叠量
        overlap = p1.radius + p2.radius - distance
        if overlap <= 0:
            return np.zeros(3)
        
        # 单位法向量
        if distance > 0:
            normal = r_vec / distance
        else:
            normal = np.array([1, 0, 0])  # 防止除以零
        
        # 相对速度
        rel_velocity = p2.velocity - p1.velocity
        
        # 法向相对速度
        vn = np.dot(rel_velocity, normal)
        
        # 赫兹接触模型参数
        R_eff = (p1.radius * p2.radius) / (p1.radius + p2.radius)
        E_eff = self.youngs_modulus / (2 * (1 - self.poissons_ratio**2))
        
        # 法向弹性力（赫兹模型）
        kn = (4/3) * E_eff * math.sqrt(R_eff * overlap)
        fn_elastic = kn * overlap
        
        # 法向阻尼力
        mass_eff = (p1.mass * p2.mass) / (p1.mass + p2.mass)
        damping_coef = -2 * math.sqrt(5/6) * self.restitution_coef * math.sqrt(mass_eff * kn)
        fn_damping = damping_coef * vn
        
        # 总法向力
        fn_total = max(fn_elastic + fn_damping, 0)  # 防止拉力
        
        # 切向力（简化的库仑摩擦）
        tangent_vel = rel_velocity - vn * normal
        if np.linalg.norm(tangent_vel) > 0:
            tangent_dir = tangent_vel / np.linalg.norm(tangent_vel)
            ft_max = self.friction_coef * fn_total
            ft_mag = min(ft_max, 0.1 * fn_total)  # 简化的切向力模型
            ft_total = -ft_mag * tangent_dir
        else:
            ft_total = np.zeros(3)
        
        # 总接触力
        contact_force = fn_total * normal + ft_total
        
        return contact_force
    
    def apply_boundary_conditions(self, particle, box_size=1.0):
        """应用边界条件"""
        boundary_force = np.zeros(3)
        boundary_stiffness = 1e6
        
        # 检查每个方向的边界
        for i in range(3):
            # 下边界
            if particle.position[i] - particle.radius < -box_size/2:
                overlap = particle.radius - (particle.position[i] + box_size/2)
                boundary_force[i] = boundary_stiffness * overlap
                # 阻尼
                boundary_force[i] -= 0.1 * particle.velocity[i]
            
            # 上边界
            if particle.position[i] + particle.radius > box_size/2:
                overlap = particle.radius - (box_size/2 - particle.position[i])
                boundary_force[i] = -boundary_stiffness * overlap
                # 阻尼
                boundary_force[i] -= 0.1 * particle.velocity[i]
        
        return boundary_force
    
    def step(self):
        """执行一个时间步长的计算"""
        # 重置所有力
        for particle in self.particles:
            particle.force = np.zeros(3)
        
        # 计算接触力
        n_particles = len(self.particles)
        for i in range(n_particles):
            for j in range(i + 1, n_particles):
                p1 = self.particles[i]
                p2 = self.particles[j]
                
                contact_force = self.calculate_contact_force(p1, p2)
                p1.force -= contact_force
                p2.force += contact_force
        
        # 应用其他力（重力、边界条件等）
        for particle in self.particles:
            if not particle.fixed:
                # 重力
                particle.force += particle.mass * self.gravity
                
                # 边界条件
                boundary_force = self.apply_boundary_conditions(particle)
                particle.force += boundary_force
        
        # 更新位置和速度（显式欧拉法）
        for particle in self.particles:
            if not particle.fixed:
                acceleration = particle.force / particle.mass
                particle.velocity += acceleration * self.dt
                particle.position += particle.velocity * self.dt
    
    def run_simulation(self, steps, visualize_every=100):
        """运行模拟"""
        print("开始模拟...")
        
        # 存储用于可视化的数据
        #positions_history = []
        
        for step in range(steps):
            self.step()
            
            if step % visualize_every == 0:
                # 记录当前位置
                current_positions = [p.position.copy() for p in self.particles]
                #positions_history.append(current_positions)
                self.visualize(current_positions, step)
                # 打印进度
                if step % (visualize_every * 10) == 0:
                    print(f"进度: {step}/{steps} 步")
        
        print("模拟完成!")
        #return positions_history
    
    def visualize(self, positions_history, time,save_path=None):
        """将每帧粒子位置与类型写入 VTK 文件（不渲染）"""

        if not positions_history:
            # 如果没有历史数据，则使用当前粒子位置作为单帧输出
            positions_history = [[p.position.copy() for p in self.particles]]

        n_particles = len(self.particles)
        if n_particles == 0:
            print("没有粒子可写入。")
            return

        # 构建每个粒子的类型（所属 cube 的索引）
        particle_types = []
        for ci, cube in enumerate(self.cubes):
            for _ in cube.particles:
                particle_types.append(ci)

        # 保证长度与粒子数匹配
        if len(particle_types) < n_particles:
            particle_types += [particle_types[-1]] * (n_particles - len(particle_types))
        elif len(particle_types) > n_particles:
            particle_types = particle_types[:n_particles]

        # 输出文件名前缀
        if save_path:
            base, ext = os.path.splitext(save_path)
            prefix = base
        else:
            prefix = "dem_particles"

        try:
            #for t, positions in enumerate(positions_history):
            if len(positions_history) != n_particles:
                    # 跳过不匹配的帧
                print(f"跳过第 {time} 帧：粒子数不匹配 ({len(positions_history)} != {n_particles})")
            

            # 创建 VTK 点和 polydata
            points = vtk.vtkPoints()
            points.SetNumberOfPoints(n_particles)
            poly = vtk.vtkPolyData()
            poly.SetPoints(points)

            # 填充点坐标
            for i, pos in enumerate(positions_history):
                points.SetPoint(i, float(pos[0]), float(pos[1]), float(pos[2]))
            points.Modified()

            # 创建整型数组表示粒子类型（所属 cube 索引）
            type_array = vtk.vtkIntArray()
            type_array.SetNumberOfComponents(1)
            type_array.SetName("Type")
            type_array.SetNumberOfTuples(n_particles)
            for i, tp in enumerate(particle_types):
                type_array.SetTuple1(i, int(tp))

            # 将类型数组附加到点数据上（标量）
            poly.GetPointData().AddArray(type_array)
            poly.GetPointData().SetScalars(type_array)

            # 写入 VTK XML PolyData (.vtp)
            filename = f"{prefix}_frame_{time:04d}.vtp"
            writer = vtk.vtkXMLPolyDataWriter()
            writer.SetFileName(filename)
            # 使用 ASCII 更易于调试/查看
            if hasattr(writer, "SetDataModeToAscii"):
                writer.SetDataModeToAscii()
            writer.SetInputData(poly)
            if writer.Write():
                print(f"Wrote {filename}")
            else:
                print(f"Failed to write {filename}")
        except Exception as e:
            print("写入 VTK 文件时出错:", e)
            return

def main():
    """主函数"""
    # 创建模拟器
    print("创建 DEM 模拟器...")
    simulator = DEMSimulator()
    
    
    # 创建多个正方体铝块
    print("创建正方体铝块...")
    
    # 正方体1 - 左侧，向右移动
    cube1 = simulator.add_cube(
        center=[-0.12, 0, 0],
        size=0.2,
        particle_radius=0.02,
        density=2700,
        velocity=[2.0, 0, 0]
    )
    
    # 正方体2 - 右侧，向左移动
    cube2 = simulator.add_cube(
        center=[0.12, 0, 0],
        size=0.2,
        particle_radius=0.02,
        density=2700,
        velocity=[-2.0, 0, 0]
    )
    
    # 正方体3 - 上方，自由落体
    '''
    cube3 = simulator.add_cube(
        center=[0, 0.3, 0.2],
        size=0.15,
        particle_radius=0.02,
        density=2700,
        velocity=[0, 0, 0]
    )
    '''
    
    print(f"总粒子数: {len(simulator.particles)}")
    
    
    # 运行模拟
    total_steps = 12000
    
    positions_history = simulator.run_simulation(
        steps=total_steps, 
        visualize_every=20
    )
    

if __name__ == "__main__":
    main()
