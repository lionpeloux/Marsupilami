using Rhino.Geometry;
using Rhino.Geometry.Intersect;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Marsupilami.Kernel
{
    public enum ConstraintTypes
    {
        Kinematic,
        Mechanical,
        Link
    }
    public enum KinematicConstraintTypes
    {
        Spherical,
        Cylindrical,
        Planar,
        PlanarOneWay,
        NurbsCrv,
        NurbsSrf,
    }
    public enum MechanicalConstraintTypes
    {
        Spring,
        PointAttractor,
        PlaneAttractor,
        SurfaceAttractor
    }
    public enum LinkConstraintTypes
    {
        SingleJoint,
        ExcentricSingleJoint,
        MultipleJoint
    }

    public abstract class DRConstraint
    {
        public ConstraintTypes constraintType;
        public abstract override string ToString(); // A UTILIER POUR AFFICHER LES INFOS D'UNE CONTRAINTE 

        // STATIC METHODES
        public static void GetConstraints(List<DRConstraint> constraints, ref DRKinematicConstraint[] kinematicConstraints, ref DRMechanicalConstraint[] mechanicalConstraints, ref DRLinkConstraint[] linkConstraints)
        {
            List<DRKinematicConstraint> listKC = new List<DRKinematicConstraint>();
            List<DRMechanicalConstraint> listMC = new List<DRMechanicalConstraint>();
            List<DRLinkConstraint> listLC = new List<DRLinkConstraint>();

            foreach (DRConstraint constraint in constraints)
            {
                switch (constraint.constraintType)
                {
                    case ConstraintTypes.Kinematic:
                        listKC.Add((DRKinematicConstraint)constraint);
                        break;
                    case ConstraintTypes.Mechanical:
                        listMC.Add((DRMechanicalConstraint)constraint);
                        break;
                    case ConstraintTypes.Link:
                        listLC.Add((DRLinkConstraint)constraint);
                        break;
                    default:
                        break;
                }
            }

            kinematicConstraints = listKC.ToArray();
            mechanicalConstraints = listMC.ToArray();
            linkConstraints = listLC.ToArray();

            Rhino.RhinoApp.WriteLine("KCST = " + kinematicConstraints.Length);
            Rhino.RhinoApp.WriteLine("MCST = " + mechanicalConstraints.Length);
            Rhino.RhinoApp.WriteLine("LCST = " + linkConstraints.Length);
        }
    }

    public abstract class DRKinematicConstraint : DRConstraint
    {
        public KinematicConstraintTypes kinematicConstraintType;
        public int ei, nj; // ADRESSE DU NOEUD

        private DRKinematicConstraint(int ei, int nj)
        {
            this.constraintType = ConstraintTypes.Kinematic;
            this.ei = ei;
            this.nj = nj;
        }
        public abstract void Enforce_V(ref Vector3d[][] Xt, ref Vector3d[][] Vt);
        public abstract void Enforce_X(ref Vector3d[][] Xt, ref Vector3d[][] Vt);

        // STATIC CONSTRUCTORS
        public static DRKinematicConstraint CreateSphericalConstraint(int ei, int nj)
        {
            DRKinematicConstraint kc = new KCSpherical(ei, nj);
            return kc;
        }
        public static DRKinematicConstraint CreateCylindricalConstraint(int ei, int nj, Vector3d u)
        {
            DRKinematicConstraint kc = new KSCylindrical(ei, nj, u);
            return kc;
        }
        public static DRKinematicConstraint CreatePlanarConstraint(int ei, int nj, Vector3d n)
        {
            DRKinematicConstraint kc = new KSPlanar(ei, nj, n);
            return kc;
        }
        public static DRKinematicConstraint CreatePlanarOneWayConstraint(int ei, int nj, Plane plane)
        {
            DRKinematicConstraint kc = new KSPlanarOneWay(ei, nj, plane);
            return kc;
        }
        public static DRKinematicConstraint CreateNurbsCrvConstraint(int ei, int nj, Curve crv)
        {
            DRKinematicConstraint kc = new KSNurbsCrv(ei, nj, crv);
            return kc;
        }
        public static DRKinematicConstraint CreateNurbsSrfConstraint(int ei, int nj, Surface srf)
        {
            DRKinematicConstraint kc = new KSNurbsSrf(ei, nj, srf);
            return kc;
        }

        // SUB-CLASS
        private class KCSpherical : DRKinematicConstraint
        {
            public KCSpherical(int ei, int nj)
                : base(ei, nj)
            {
                this.kinematicConstraintType = KinematicConstraintTypes.Spherical;

            }
            public override string ToString()
            {
                return "KinematicConstraint[Spherical] = [" + ei + "," + nj + "]";
            }
            public override void Enforce_V(ref Vector3d[][] Xt, ref Vector3d[][] Vt)
            {
                Vt[ei][nj].X = 0.0; Vt[ei][nj].Y = 0.0; Vt[ei][nj].Z = 0.0;
            }
            public override void Enforce_X(ref Vector3d[][] Xt, ref Vector3d[][] Vt)
            {
            }
        }
        private class KSCylindrical : DRKinematicConstraint
        {
            protected Vector3d u; // VECTEUR DIRECTEUR DU SLIDER 

            public KSCylindrical(int ei, int nj, Vector3d u)
                : base(ei, nj)
            {
                this.kinematicConstraintType = KinematicConstraintTypes.Cylindrical;
                this.u = u;
                this.u.Unitize();
            }
            public override string ToString()
            {
                return "KinematicConstraint[Cylindrical] = [" + ei + "," + nj + "]";
            }
            public override void Enforce_V(ref Vector3d[][] Xt, ref Vector3d[][] Vt)
            {
                Vt[ei][nj] = Vector3d.Multiply(Vt[ei][nj], this.u) * this.u;
            }
            public override void Enforce_X(ref Vector3d[][] Xt, ref Vector3d[][] Vt)
            {
            }
        }
        private class KSPlanar : DRKinematicConstraint
        {
            protected Vector3d n; // VECTEUR NORMAL AU PLAN DU SLIDER 

            public KSPlanar(int ei, int nj, Vector3d n)
                : base(ei, nj)
            {
                this.kinematicConstraintType = KinematicConstraintTypes.Planar;
                this.n = n;
                this.n.Unitize();
            }
            public override string ToString()
            {
                return "KinematicConstraint[Planar] = [" + ei + "," + nj + "]";
            }
            public override void Enforce_V(ref Vector3d[][] Xt, ref Vector3d[][] Vt)
            {
                Vt[ei][nj] = Vt[ei][nj] - Vector3d.Multiply(Vt[ei][nj], this.n) * this.n;
            }
            public override void Enforce_X(ref Vector3d[][] Xt, ref Vector3d[][] Vt)
            {
            }
        }
        private class KSPlanarOneWay : DRKinematicConstraint
        {
            protected Plane plane;        // LE PLAN EST ORIENTE DANS LE SENS DU DOF
            protected Point3d Xprev;      // POSITION PRECEDENTE POUR DETECTER LE CONTACT
            protected Line trajectory;    // TRAJECTOIRE ENTRE X(t) et X(t+dt)
            protected double param;
            protected double d = 1e-5;

            public KSPlanarOneWay(int ei, int nj, Plane plane)
                : base(ei, nj)
            {
                this.kinematicConstraintType = KinematicConstraintTypes.Planar;
                this.plane = plane;
            }
            public override string ToString()
            {
                return "KinematicConstraint[PlanarOneWay] = [" + ei + "," + nj + "]";
            }
            public override void Enforce_V(ref Vector3d[][] Xt, ref Vector3d[][] Vt)
            {
                Xprev = new Point3d(Xt[ei][nj]);

                if (plane.DistanceTo(Xprev) < d)
                {
                    if (Vector3d.Multiply(Vt[ei][nj], plane.Normal) < 0)
                    {
                        Vt[ei][nj] = Vt[ei][nj] - Vector3d.Multiply(Vt[ei][nj], plane.Normal) * plane.Normal;
                    }
                }
            }
            public override void Enforce_X(ref Vector3d[][] Xt, ref Vector3d[][] Vt)
            {
                trajectory = new Line(Xprev, new Point3d(Xt[ei][nj]));
                Intersection.LinePlane(trajectory, plane, out param);

                if (0 <= param && param <= 1) // la trajectoire intersecte le plan
                {
                    Xt[ei][nj] = (Vector3d)trajectory.PointAt(param);
                }
            }
        }
        private class KSNurbsCrv : DRKinematicConstraint
        {
            protected Curve crv;  // VECTEUR NORMAL AU PLAN DU PLANAR SLIDER
            protected Vector3d t;
            protected double param;

            public KSNurbsCrv(int ei, int nj, Curve crv)
                : base(ei, nj)
            {
                this.kinematicConstraintType = KinematicConstraintTypes.NurbsCrv;
                this.crv = crv;
            }
            public override string ToString()
            {
                return "KinematicConstraint[NurbsCrv] = [" + ei + "," + nj + "]";
            }
            public override void Enforce_V(ref Vector3d[][] Xt, ref Vector3d[][] Vt)
            {
                crv.ClosestPoint((Point3d)Xt[ei][nj], out param);
                t = crv.TangentAt(param);
                Vt[ei][nj] = Vector3d.Multiply(Vt[ei][nj], t) * t;
            }
            public override void Enforce_X(ref Vector3d[][] Xt, ref Vector3d[][] Vt)
            {
                crv.ClosestPoint((Point3d)Xt[ei][nj], out param);
                Xt[ei][nj] = (Vector3d)crv.PointAt(param);
            }
        }
        private class KSNurbsSrf : DRKinematicConstraint
        {
            protected Surface srf;    // SURFACE DU SLIDER
            protected double param_u; // PARAMETRE U
            protected double param_v; // PARAMETRE V
            protected Vector3d n;     // NORMALE LOCAL

            public KSNurbsSrf(int ei, int nj, Surface srf)
                : base(ei, nj)
            {
                this.kinematicConstraintType = KinematicConstraintTypes.NurbsCrv;
                this.srf = srf;
            }
            public override string ToString()
            {
                return "KinematicConstraint[NrbsSrf] = [" + ei + "," + nj + "]";
            }
            public override void Enforce_V(ref Vector3d[][] Xt, ref Vector3d[][] Vt)
            {
                srf.ClosestPoint((Point3d)Xt[ei][nj], out param_u, out param_v);
                n = srf.NormalAt(param_u, param_v);
                Vt[ei][nj] = Vt[ei][nj] - Vector3d.Multiply(Vt[ei][nj], n) * n;
            }
            public override void Enforce_X(ref Vector3d[][] Xt, ref Vector3d[][] Vt)
            {
                srf.ClosestPoint((Point3d)Xt[ei][nj], out param_u, out param_v);
                Xt[ei][nj] = (Vector3d)srf.PointAt(param_u, param_v);
            }
        }
    }
    public abstract class DRMechanicalConstraint : DRConstraint
    {
        public MechanicalConstraintTypes mechanicalConstraintType;
        protected int ei, nj;   // ADRESSE DU NOEUD

        private DRMechanicalConstraint(int ei, int nj)
        {
            this.constraintType = ConstraintTypes.Mechanical;
            this.ei = ei;
            this.nj = nj;
        }
        public abstract void Update_LM(ref double[][] LMt, ref Vector3d[][] Rt);
        public abstract void Update_R(ref Vector3d[][] Xt, ref Vector3d[][] Rt);

        // STATIC CONSTRUCTORS
        public static DRMechanicalConstraint CreateSpringConstraint(int ei, int nj, Vector3d Xpin, double k)
        {
            DRMechanicalConstraint mc = new MCSpring(ei, nj, Xpin, k, k, k);
            return mc;
        }
        public static DRMechanicalConstraint CreateSpringConstraint(int ei, int nj, Vector3d Xpin, double kx, double ky, double kz)
        {
            DRMechanicalConstraint mc = new MCSpring(ei, nj, Xpin, kx, ky, kz);
            return mc;
        }
        public static DRMechanicalConstraint CreatePointAttractorConstraint(int ei, int nj, Vector3d Xpin, double k)
        {
            DRMechanicalConstraint mc = new MCPointAttractor(ei, nj, Xpin, k);
            return mc;
        }
        public static DRMechanicalConstraint CreatePlaneAttractorConstraint(int ei, int nj, Plane plane, double k)
        {
            DRMechanicalConstraint mc = new MCPlaneAttractor(ei, nj, plane, k);
            return mc;
        }
        public static DRMechanicalConstraint CreateSurfaceAttractorConstraint(int ei, int nj, Surface srf, double k)
        {
            DRMechanicalConstraint mc = new MCSurfaceAttractor(ei, nj, srf, k);
            return mc;
        }

        // SUB-CLASS
        private class MCSpring : DRMechanicalConstraint
        {
            protected Vector3d Xpin;      // POSITION D'ACCROCHE
            protected double kx, ky, kz;  // RAIDEUR DU RESSORT
            protected Vector3d u;         // VECTEUR DIRECTEUR
            protected Vector3d F;         // EFFORT DE RAPPEL SUR LE NOEUD [ei,nj]

            public MCSpring(int ei, int nj, Vector3d Xpin, double kx, double ky, double kz)
                : base(ei, nj)
            {
                this.mechanicalConstraintType = MechanicalConstraintTypes.PointAttractor;
                this.Xpin = Xpin;
                this.kx = kx;
                this.ky = ky;
                this.kz = kz;
            }
            public override string ToString()
            {
                return "MechanicalConstraint[Spring] = [" + ei + "," + nj + "]" + " | " + string.Format("kx = {0:F3}, ky = {1:F3}, kz = {2:F3}", kx, ky, kz);
            }
            public override void Update_LM(ref double[][] LMt, ref Vector3d[][] Rt)
            {
            }
            public override void Update_R(ref Vector3d[][] Xt, ref Vector3d[][] Rt)
            {
                u = Xpin - Xt[ei][nj];
                F = new Vector3d(kx * u.X, ky * u.Y, kz * u.Z);
                Rt[ei][nj] += F;
            }
        }
        private class MCPointAttractor : DRMechanicalConstraint
        {
            protected Vector3d Xpin;  // POSITION D'ACCROCHE    
            protected double k;       // RAIDEUR DU RESSORT
            protected Vector3d F;     // EFFORT DE RAPPEL SUR LE NOEUD [ei,nj]

            public MCPointAttractor(int ei, int nj, Vector3d Xpin, double k)
                : base(ei, nj)
            {
                this.mechanicalConstraintType = MechanicalConstraintTypes.PointAttractor;
                this.Xpin = Xpin;
                this.k = k;
            }
            public override string ToString()
            {
                return "MechanicalConstraint[PointAttractor] = [" + ei + "," + nj + "]" + " | " + string.Format("kx = {0:F3}, ky = {1:F3}, kz = {2:F3}", k, k, k);
            }
            public override void Update_LM(ref double[][] LMt, ref Vector3d[][] Rt)
            {
            }
            public override void Update_R(ref Vector3d[][] Xt, ref Vector3d[][] Rt)
            {
                F = k * (Xpin - Xt[ei][nj]);
                Rt[ei][nj] += F;
            }
        }
        private class MCPlaneAttractor : DRMechanicalConstraint
        {
            protected Plane plane;    // ATTRACTEUR    
            protected double k;       // RAIDEUR DU RESSORT
            protected Vector3d Xproj; // PROJECTION
            protected Vector3d F;     // EFFORT DE RAPPEL SUR LE NOEUD [ei,nj]

            public MCPlaneAttractor(int ei, int nj, Plane plane, double k)
                : base(ei, nj)
            {
                this.mechanicalConstraintType = MechanicalConstraintTypes.PlaneAttractor;
                this.plane = plane;
                this.k = k;
            }
            public override string ToString()
            {
                return "MechanicalConstraint[PlaneAttractor] = [" + ei + "," + nj + "]" + " | " + string.Format("k = {0:F3}", k);
            }
            public override void Update_LM(ref double[][] LMt, ref Vector3d[][] Rt)
            {
            }
            public override void Update_R(ref Vector3d[][] Xt, ref Vector3d[][] Rt)
            {
                Xproj = (Vector3d)plane.ClosestPoint((Point3d)Xt[ei][nj]);
                F = k * (Xproj - Xt[ei][nj]);
                Rt[ei][nj] += F;
            }
        }
        private class MCSurfaceAttractor : DRMechanicalConstraint
        {
            protected Surface srf;    // ATTRACTEUR    
            protected double k;       // RAIDEUR DU RESSORT
            protected double param_u; // PARAMETRE U DU PROJETE
            protected double param_v; // PARAMETRE V DU PROJETE
            protected Vector3d Xproj; // PROJECTION
            protected Vector3d F;     // EFFORT DE RAPPEL SUR LE NOEUD [ei,nj]

            public MCSurfaceAttractor(int ei, int nj, Surface srf, double k)
                : base(ei, nj)
            {
                this.mechanicalConstraintType = MechanicalConstraintTypes.SurfaceAttractor;
                this.srf = srf;
                this.k = k;
            }
            public override string ToString()
            {
                return "MechanicalConstraint[PlaneAttractor] = [" + ei + "," + nj + "]" + " | " + string.Format("k = {0:F3}", k);
            }
            public override void Update_LM(ref double[][] LMt, ref Vector3d[][] Rt)
            {
            }
            public override void Update_R(ref Vector3d[][] Xt, ref Vector3d[][] Rt)
            {
                if (srf.ClosestPoint((Point3d)Xt[ei][nj], out param_u, out param_v))
                {
                    Xproj = (Vector3d)srf.PointAt(param_u, param_v);
                    F = k * (Xproj - Xt[ei][nj]);
                    Rt[ei][nj] += F;
                }
            }
        }
    }
    public abstract class DRLinkConstraint : DRConstraint
    {
        public LinkConstraintTypes linkConstraintType;

        protected int count;        // NOMBRE DE NOEUDS DU LINK
        protected Vector3d R;       // EFFORT RESULTANT GLOBAL
        protected double lm;        // LUMPED MASS RESULTANTE GLOBALE
        protected Vector3d Xbar;    // POSITION BARYCENTRIQUE DU LINK

        private DRLinkConstraint()
        {
            this.constraintType = ConstraintTypes.Link;
        }
        public abstract void Transfer_LM(ref double[][] LMt);
        public abstract void Transfer_R(ref Vector3d[][] Rt);
        public abstract void Enforce_X(ref Vector3d[][] Xt);

        // STATIC CONSTRUCTORS
        public static DRLinkConstraint CreateSingleJointConstraint(int ei1, int nj1, int ei2, int nj2)
        {
            DRLinkConstraint lc = new LCSingleJoint(ei1, nj1, ei2, nj2);
            return lc;
        }
        public static DRLinkConstraint CreateExcentricSingleJointConstraint(int ei1, int nj1, int ei2, int nj2, double e)
        {
            DRLinkConstraint lc = new LCExcentricSingleJoint(ei1, nj1, ei2, nj2, e);
            return lc;
        }
        public static DRLinkConstraint CreateMultipleJointConstraint(int[] ei, int[] nj)
        {
            DRLinkConstraint lc = new LCMultipleJoint(ei, nj);
            return lc;
        }

        // SUB-CLASS
        private class LCSingleJoint : DRLinkConstraint
        {
            protected int ei1, nj1;        // ADRESSE DU NOEUD 1
            protected int ei2, nj2;        // ADRESSE DU NOEUD 2

            public LCSingleJoint(int ei1, int nj1, int ei2, int nj2)
                : base()
            {
                this.linkConstraintType = LinkConstraintTypes.SingleJoint;
                this.count = 2;
                this.ei1 = ei1;
                this.nj1 = nj1;
                this.ei2 = ei2;
                this.nj2 = nj2;
            }
            public override string ToString()
            {
                return "LinkConstraint[SingleJoint] = " + "[" + ei1 + "," + nj1 + "]" + " | " + "[" + ei2 + "," + nj2 + "]" + "]";
            }

            public override void Transfer_LM(ref double[][] LMt)
            {
                lm = LMt[ei1][nj1] + LMt[ei2][nj2];
                LMt[ei1][nj1] = lm;
                LMt[ei2][nj2] = lm;
            }
            public override void Transfer_R(ref Vector3d[][] Rt)
            {
                R = Rt[ei1][nj1] + Rt[ei2][nj2];
                Rt[ei1][nj1] = R;
                Rt[ei2][nj2] = R;
            }
            public override void Enforce_X(ref Vector3d[][] Xt)
            {
                Xbar = (Xt[ei1][nj1] + Xt[ei2][nj2]) / 2;

                Xt[ei1][nj1] = Xbar;
                Xt[ei2][nj2] = Xbar;
            }
        }
        private class LCExcentricSingleJoint : LCSingleJoint
        {
            protected double e;           // EXCENTRICITE
            protected Vector3d X1, X2;    // POSITIONS EXCENTREE
            protected Vector3d t1, t2;    // VECTEURS DIRECTEURS LOCAUX DES ELEMENTS
            protected Vector3d n;         // VECTEUR DIRECTEUR DE LA LIAISON

            public LCExcentricSingleJoint(int ei1, int nj1, int ei2, int nj2, double e)
                : base(ei1, nj1, ei2, nj2)
            {
                this.linkConstraintType = LinkConstraintTypes.ExcentricSingleJoint;
                this.count = 2;
                this.ei1 = ei1;
                this.nj1 = nj1;
                this.ei2 = ei2;
                this.nj2 = nj2;
                this.e = e;
            }
            public override string ToString()
            {
                return "LinkConstraint[ExcentricSingleJoint] = " + "[" + ei1 + "," + nj1 + "]" + " | " + "[" + ei2 + "," + nj2 + "]" + "]";
            }

            public override void Transfer_LM(ref double[][] LMt)
            {
                lm = LMt[ei1][nj1] + LMt[ei2][nj2];
                LMt[ei1][nj1] = lm;
                LMt[ei2][nj2] = lm;
            }
            public override void Transfer_R(ref Vector3d[][] Rt)
            {
                R = Rt[ei1][nj1] + Rt[ei2][nj2];
                Rt[ei1][nj1] = R;
                Rt[ei2][nj2] = R;
            }
            public override void Enforce_X(ref Vector3d[][] Xt)
            {
                Xbar = (Xt[ei1][nj1] + Xt[ei2][nj2]) / 2;

                t1 = (Xt[ei1][nj1 - 1] - Xt[ei1][nj1 + 1]);
                t2 = (Xt[ei2][nj2 - 1] - Xt[ei2][nj2 + 1]);
                n = Vector3d.CrossProduct(t1, t2);
                n.Unitize();

                X1 = Xbar + e / 2 * n;
                X2 = Xbar - e / 2 * n;

                Xt[ei1][nj1] = X1;
                Xt[ei2][nj2] = X2;
            }
        }
        private class LCMultipleJoint : DRLinkConstraint
        {
            protected int[] ei, nj;       // ADRESSE DES NOEUDS

            public LCMultipleJoint(int[] ei, int[] nj)
                : base()
            {
                this.linkConstraintType = LinkConstraintTypes.MultipleJoint;
                this.count = ei.Length;
                this.ei = ei;
                this.nj = nj;
            }
            public override string ToString()
            {
                string str = "LinkConstraint[MultipleJoint] = ";
                for (int k = 0; k < count - 1; k++) { str += "[" + ei[k] + "," + nj[k] + "]" + " | "; }
                str += "[" + ei[count - 1] + "," + nj[count - 1] + "]";
                return str;
            }

            public override void Transfer_LM(ref double[][] LMt)
            {
                lm = 0;
                for (int k = 0; k < count; k++)
                {
                    lm += LMt[ei[k]][nj[k]];
                }

                for (int k = 0; k < count; k++)
                {
                    LMt[ei[k]][nj[k]] = lm;
                }
            }
            public override void Transfer_R(ref Vector3d[][] Rt)
            {
                R = new Vector3d(0, 0, 0);
                for (int k = 0; k < count; k++)
                {
                    R += Rt[ei[k]][nj[k]];
                }

                for (int k = 0; k < count; k++)
                {
                    Rt[ei[k]][nj[k]] = R;
                }
            }
            public override void Enforce_X(ref Vector3d[][] Xt)
            {
                Xbar = new Vector3d(0, 0, 0);

                for (int k = 0; k < count; k++)
                {
                    Xbar += Xt[ei[k]][nj[k]];
                }

                Xbar = Xbar / count;

                for (int k = 0; k < count; k++)
                {
                    Xt[ei[k]][nj[k]] = Xbar;
                }
            }
        }
    }
}
