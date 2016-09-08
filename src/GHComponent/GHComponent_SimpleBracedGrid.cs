using Grasshopper.Kernel;
using Marsupilami.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace Marsupilami.GHComponent
{
    public class GHComponent_SimpleBracedGrid : GH_Component
    {
        
        #region DYNAMIC OUTPUTS ENGINE
        // specific code to refresh outputs
        // FLAGS
        public enum LoopFlag
        {
            running,    // process is running
            aborted,    // process has been manually aborted by the user pressing the esc key
            stoped,     // process has been stoped before iteration_max was reached
            ended,      // process has reached iteration_max
        }
        private int iteration_current = 1;
        private int iteration_groupby;      // group of iterations between output refresh
        private int iteration_max;          // total iteration to run
        private bool loop_reset = true;
        private bool loop_reset_cache = true;
        private LoopFlag flag = LoopFlag.ended;
        private GH_Document ghDocument;
        private System.Windows.Forms.Timer trampoline = new System.Windows.Forms.Timer(); // needed to restart solution for intermediate outputs refresh
        // METHODS
        private void documentSolutionEnd(object sender, GH_SolutionEventArgs e)
        {
            ghDocument.SolutionEnd -= documentSolutionEnd;

            // The trampoline gives time to stop the infinitly recursive chain from canvas,
            // as well as mantains the call stack fixed in length
            trampoline.Interval = 1;
            trampoline.Tick += trampolineTick;
            trampoline.Start();
        }
        private void trampolineTick(object sender, EventArgs e)
        {
            trampoline.Tick -= trampolineTick;
            trampoline.Stop();
            this.ExpireSolution(true);
        }
        #endregion

        // FIELDS
        private int index;
        int nH, nV, n1, n2;
        private List<Curve> crvs_H, crvs_V, crvs_1, crvs_2;
        private Polyline[] polylines_H, polylines_V, polylines_1, polylines_2;
        private List<Point3d> pts;
        private double Eg, Sg, Ig, Eb, Sb, Ib, fz, dij;
       


        public DRElement[] elements;
        public List<DRConstraint> constraints;
        public DRRelax solver;
        List<double> Ec;
        public double dt;
        public double Eclim;

        // CONSTRUCTOR
        public GHComponent_SimpleBracedGrid()
            : base("Simple Braced Grid", "Simple Braced Grid", "Relaxation of a braced grid on hinged supports", "Marsupilami", "Example")
        {
            //if (DateTime.Now > MarsupilamiInfo.ExpirationDate) { throw new Exception("Enable to launch the gridshell library. Please contact the author."); }

            crvs_H = new List<Curve>();
            crvs_V = new List<Curve>();
            crvs_1 = new List<Curve>();
            crvs_2 = new List<Curve>();
            pts = new List<Point3d>();

            Ec = new List<double>();
            dt = 1;
            Eclim = 1e-10;

            
        }
        public override Guid ComponentGuid
        {
            get { return new Guid("{a9d1618d-39db-4672-b1b9-2e66dbedf194}"); }
        }
        public override GH_Exposure Exposure
        {
            get { return GH_Exposure.primary; }
        }

        // PARAMETERS
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Polylines H", "H", "List of polylines in the horizontal direction", GH_ParamAccess.list);
            pManager.AddCurveParameter("Polylines V", "V", "List of polylines in the vertical direction", GH_ParamAccess.list);
            pManager.AddCurveParameter("Polylines D1", "D1", "List of polylines in the first bracing direction", GH_ParamAccess.list);
            pManager.AddCurveParameter("Polylines D2", "D2", "List of polylines in the second bracing direction", GH_ParamAccess.list);

            pManager.AddPointParameter("Hinged Point", "pts", "List of hinged supports locations", GH_ParamAccess.list);

            pManager.AddNumberParameter("Young's Modulus", "Eg", "Young's modulus of the beam's material (Pa)", GH_ParamAccess.item, 25e9);
            pManager.AddNumberParameter("Cross Section Area", "Sg", "Area of the beam's cross section (m2)", GH_ParamAccess.item, 4.2e-4);
            pManager.AddNumberParameter("Area Moment of Inertia", "Ig", "Area moment of inertia of the beam's cross section (m4)", GH_ParamAccess.item, 7.9e-8);
            pManager.AddNumberParameter("Young's Modulus", "Eb", "Young's modulus of the bracing's material (Pa)", GH_ParamAccess.item, 25e9);
            pManager.AddNumberParameter("Cross Section Area", "Sb", "Area of the bracing's cross section (m2)", GH_ParamAccess.item, 4.2e-4);
            pManager.AddNumberParameter("Area Moment of Inertia", "Ib", "Area moment of inertia of the bracing's cross section (m4)", GH_ParamAccess.item, 7.9e-8);

            pManager.AddNumberParameter("Vertical Force", "fz", "Vertical force to apply per node", GH_ParamAccess.item, 0);

            pManager.AddIntegerParameter("iteration_max", "N_max", "Total number of iterations to run", GH_ParamAccess.item, 1000);
            pManager.AddIntegerParameter("iterations_groupby", "N_refresh", "Group of iterations between output refresh", GH_ParamAccess.item, 1000);
            pManager.AddBooleanParameter("reset", "reset", "Reset the engine", GH_ParamAccess.item);

            pManager[0].Optional = false;
            pManager[1].Optional = false;
            pManager[2].Optional = false;
            pManager[3].Optional = false;
            pManager[4].Optional = false;
            pManager[5].Optional = false;
            pManager[6].Optional = false;
            pManager[7].Optional = false;
            pManager[8].Optional = false;
            pManager[9].Optional = false;
            pManager[10].Optional = false;
            pManager[11].Optional = false;
            pManager[12].Optional = false;
            pManager[13].Optional = false;
            pManager[14].Optional = false;
        }
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.Register_StringParam("info", "info", "");
            pManager.Register_PointParam("Xi", "Xi", "");
            pManager.Register_VectorParam("Xt", "Xt", "");
            pManager.Register_VectorParam("Vt", "Vt", "");
            pManager.Register_DoubleParam("LMt", "LMt", "");
            pManager.Register_VectorParam("Rt", "Rt", "");
            pManager.Register_VectorParam("Fext", "Fext", "");
            pManager.Register_DoubleParam("Ect", "Ect", "");
            pManager.Register_DoubleParam("Nt", "Nt", "");
            pManager.Register_DoubleParam("Tt", "Tt", "");
            pManager.Register_DoubleParam("Mt", "Mt", "");
        }

        // SOLVER
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            ghDocument = this.OnPingDocument();
        
            pts.Clear();
          
            if (!DA.GetDataList(4, pts)) { return; }

            if (!DA.GetData(5, ref Eg)) { return; }
            if (!DA.GetData(6, ref Sg)) { return; }
            if (!DA.GetData(7, ref Ig)) { return; }
            if (!DA.GetData(8, ref Eb)) { return; }
            if (!DA.GetData(9, ref Sb)) { return; }
            if (!DA.GetData(10, ref Ib)) { return; }

            if (!DA.GetData(11, ref fz)) { return; }

            if (!DA.GetData(12, ref iteration_max)) { return; }
            if (!DA.GetData(13, ref iteration_groupby)) { return; }
            if (!DA.GetData(14, ref loop_reset)) { return; }
               

            if (solver == null || loop_reset != loop_reset_cache) // Premier Calcul
            {
                // EXTRACT POLYLINES FROM INPUT CURVES 
                crvs_H.Clear();
                crvs_V.Clear();
                crvs_1.Clear();
                crvs_2.Clear();
                if (!DA.GetDataList(0, crvs_H)) { return; }
                if (!DA.GetDataList(1, crvs_V)) { return; }
                if (!DA.GetDataList(2, crvs_1)) { return; }
                if (!DA.GetDataList(3, crvs_2)) { return; }

                nH = crvs_H.Count;
                polylines_H = new Polyline[nH];
                for (int i = 0; i < nH; i++)
                {
                    crvs_H[i].TryGetPolyline(out polylines_H[i]);
                }
                nV = crvs_V.Count;
                polylines_V = new Polyline[nV];
                for (int j = 0; j < nV; j++)
                {
                    crvs_V[j].TryGetPolyline(out polylines_V[j]);
                }
                n1 = crvs_1.Count;
                polylines_1 = new Polyline[n1];
                for (int i = 0; i < n1; i++)
                {
                    crvs_1[i].TryGetPolyline(out polylines_1[i]);
                }
                n2 = crvs_2.Count;
                polylines_2 = new Polyline[n2];
                for (int j = 0; j < n2; j++)
                {
                    crvs_2[j].TryGetPolyline(out polylines_2[j]);
                }

                // INIT
                loop_reset_cache = loop_reset;
                Init();

                // FIRST ITERATION
                ghDocument.SolutionEnd += documentSolutionEnd;
                Run();
            }
            else
            {
                if (flag == LoopFlag.running)
                {
                    ghDocument.SolutionEnd += documentSolutionEnd;

                    // GROUP OF ITERATIONS
                    while (flag == LoopFlag.running)
                    {
                        // SINGLE ITERATION
                        iteration_current++;
                        Run();

                        // LOOPSTATE : STOP
                        if (flag == LoopFlag.stoped)
                        {
                            ghDocument.SolutionEnd -= documentSolutionEnd;
                            break;
                        }

                        // LOOPSTATE : END
                        if (iteration_current == iteration_max)
                        {
                            ghDocument.SolutionEnd -= documentSolutionEnd;
                            flag = LoopFlag.ended;
                            break;
                        }

                        // LOOP STATE : ABORT
                        if (GH_Document.IsEscapeKeyDown())
                        {
                            ghDocument.SolutionEnd -= documentSolutionEnd;
                            flag = LoopFlag.aborted;
                            break;
                        }

                        // REFRESH DETECTION
                        if (iteration_current % iteration_groupby == 0) { OnRefresh(); break; }
                    }
                }
            }

            AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, flag.ToString());

            DA.SetDataList(0, solver.Info());
            DA.SetDataTree(1, solver.GHTree_Xi());
            DA.SetDataTree(2, solver.GHTree_Xt());
            DA.SetDataTree(3, solver.GHTree_Vt());
            DA.SetDataTree(4, solver.GHTree_LMt());
            DA.SetDataTree(5, solver.GHTree_Rt());
            DA.SetDataTree(6, solver.GHTree_Fext());
            DA.SetDataList(7, Ec);
            DA.SetDataTree(8, solver.GHTree_N());
            DA.SetDataTree(9, solver.GHTree_T());
            DA.SetDataTree(10, solver.GHTree_M());
        }

        private void Init()
        {
            // reset flags, do first iteration, refresh outputs and restart iterative process
            iteration_current = 1;
            flag = LoopFlag.running;

            elements = new DRElement[nH + nV + n1 + n2];
            constraints = new List<DRConstraint>();
            
            // ELEMENT H
            for (int i = 0; i < nH; i++)
            {
                int n = polylines_H[i].Count;
                var Xi = new Vector3d[n];
                var L0 = new double[n-1];
                var Fext = new Vector3d[n];

                //var Xi = new Vector3d[2*n-1];
                //var L0 = new double[2*n-2];
                //var Fext = new Vector3d[2*n-1];

                Xi[0] = (Vector3d)polylines_H[i][0];
                Fext[0] = new Vector3d(0,0,fz);
                for (int k = 1; k < n; k++)
                {
                    Xi[k] = (Vector3d)polylines_H[i][k];
                    Fext[k] = new Vector3d(0,0,fz);
                    L0[k-1] = (Xi[k]- Xi[k-1]).Length;
                }
                //var xxi = new Vector3d(0,0,0);
                //for (int k = 1; k < n; k++)
                //{
                //    xxi = (Vector3d)polylines_H[i][k];
                //    Xi[2*k-1] = (xxi + Xi[2*k-2])/2;
                //    Xi[2*k] = xxi;
                //    Fext[2*k-1] = new Vector3d(0, 0, fz);
                //    Fext[2*k] = new Vector3d(0, 0, fz);
                //    L0[2*k-2] = (Xi[2*k-1] - Xi[2*k-2]).Length;
                //    L0[2*k-1] = (Xi[2*k] - Xi[2*k-1]).Length;
                //}
                elements[i] = DRBeamElement.CreateSimpleBeamElement(i, Xi, L0, Fext, dt);
                elements[i].E = Eg;
                elements[i].S = Sg;
                elements[i].I = Ig;
            }

            // ELEMENT V
            for (int j = nH; j < nH+nV; j++)
            {
                int n = polylines_V[j-nH].Count;
                var Xi = new Vector3d[n];
                var L0 = new double[n-1];
                var Fext = new Vector3d[n];
                
                Xi[0] = (Vector3d)polylines_V[j - nH][0];
                Fext[0] = new Vector3d(0, 0, 0);
                for (int k = 1; k < n; k++)
                {
                    Xi[k] = (Vector3d)polylines_V[j - nH][k];
                    Fext[k] = new Vector3d(0, 0, 0);
                    L0[k-1] = (Xi[k] - Xi[k - 1]).Length;
                }

                elements[j] = DRBeamElement.CreateSimpleBeamElement(j, Xi, L0, Fext, dt);
                elements[j].E = Eg;
                elements[j].S = Sg;
                elements[j].I = Ig;
            }

            // ELEMENT D1
            for (int i = nH + nV; i < n1 + nH + nV; i++)
            {
                int n = polylines_1[i-nH-nV].Count;
                var Xi = new Vector3d[n];
                var L0 = new double[n - 1];
                var Fext = new Vector3d[n];

                Xi[0] = (Vector3d)polylines_1[i - nH - nV][0];
                Fext[0] = new Vector3d(0, 0, fz);
                for (int k = 1; k < n; k++)
                {
                    Xi[k] = (Vector3d)polylines_1[i - nH - nV][k];
                    Fext[k] = new Vector3d(0, 0, fz);
                    L0[k - 1] = (Xi[k] - Xi[k - 1]).Length;
                }
                
                elements[i] = DRBeamElement.CreateSimpleBeamElement(i, Xi, L0, Fext, dt);
                elements[i].E = Eb;
                elements[i].S = Sb;
                elements[i].I = Ib;
            }

            //// ELEMENT D2
            for (int i = nH+nV+n1; i < n2 + n1 + nH + nV; i++)
            {
                int n = polylines_2[i - nH - nV - n1].Count;
                var Xi = new Vector3d[n ];
                var L0 = new double[n - 1];
                var Fext = new Vector3d[n];

                Xi[0] = (Vector3d)polylines_2[i - nH - nV - n1][0];
                Fext[0] = new Vector3d(0, 0, fz);
                for (int k = 1; k < n; k++)
                {
                    Xi[k] = (Vector3d)polylines_2[i - nH - nV - n1][k];
                    Fext[k] = new Vector3d(0, 0, fz);
                    L0[k - 1] = (Xi[k] - Xi[k - 1]).Length;
                }

                elements[i] = DRBeamElement.CreateSimpleBeamElement(i, Xi, L0, Fext, dt);
                elements[i].E = Eb;
                elements[i].S = Sb;
                elements[i].I = Ib;
            }

            // LINKS
            int[,] link = new int[4, 2];
            for (int i = 0; i < nH; i++)
            {
                for (int ki = 0; ki < polylines_H[i].Count; ki++)
                {
                    link[0, 0] = i;
                    link[0, 1] = ki;
                    for (int j = 0; j < nV; j++)
                    {
                        for (int kj = 0; kj < polylines_V[j].Count; kj++)
                        {
                            if (polylines_H[i][ki].DistanceTo(polylines_V[j][kj]) < 0.001)
                            {
                                link[1, 0] = j + nH;
                                link[1, 1] = kj;
                                dij = 0;
                                for (int u1 = 0; u1 < n1; u1++)
                                {
                                    for (int k1 = 0; k1 < polylines_1[u1].Count; k1++)
                                    {
                                        if (polylines_H[i][ki].DistanceTo(polylines_1[u1][k1]) < 0.001)
                                        {
                                            dij = 1;
                                            for (int u2 = 0; u2 < n2; u2++)
                                            {
                                                for (int k2 = 0; k2 < polylines_2[u2].Count; k2++)
                                                {
                                                    if (polylines_H[i][ki].DistanceTo(polylines_2[u2][k2]) < 0.001)
                                                    {
                                                        link[2, 0] = u1 + nH + nV;
                                                        link[2, 1] = k1;
                                                        link[3, 0] = u2 + nH + nV + n1;
                                                        link[3, 1] = k2;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                int[] ei = new int[4] { link[0, 0], link[1, 0], link[2, 0], link[3, 0] };
                                int[] nj = new int[4] { link[0, 1], link[1, 1], link[2, 1], link[3, 1] };
                                if (dij == 0)
                                {
                                    constraints.Add(DRLinkConstraint.CreateSingleJointConstraint(i, ki, j + nH, kj));
                                }
                                else
                                {
                                    constraints.Add(DRLinkConstraint.CreateMultipleJointConstraint(ei, nj));
                                }
                            }
                        }

                    }
                }
            }

            // HINGED SUPPORTS
            for (int i = 0; i < nH; i++)
            {
                for (int ki = 0; ki < polylines_H[i].Count; ki++)
                {
                    for (int j = 0; j < pts.Count; j++)
                    {
                        if (polylines_H[i][ki].DistanceTo(pts[j]) < 0.001)
                        {
                            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i, ki));
                            //constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i, 2 * ki));
                        }
                    }
                }
            }
            for (int i = 0; i < nV; i++)
            {
                for (int ki = 0; ki < polylines_V[i].Count; ki++)
                {
                    for (int j = 0; j < pts.Count; j++)
                    {
                        if (polylines_V[i][ki].DistanceTo(pts[j]) < 0.001)
                        {
                            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i + nH, ki));
                            //constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i + nH, 2 * ki));
                        }
                    }
                }
            }
            for (int i = 0; i < n1; i++)
            {
                for (int ki = 0; ki < polylines_1[i].Count; ki++)
                {
                    for (int j = 0; j < pts.Count; j++)
                    {
                        if (polylines_1[i][ki].DistanceTo(pts[j]) < 0.001)
                        {
                            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i + nH + nV, ki));
                            //constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i + nH, 2 * ki));
                        }
                    }
                }
            }
            for (int i = 0; i < n2; i++)
            {
                for (int ki = 0; ki < polylines_2[i].Count; ki++)
                {
                    for (int j = 0; j < pts.Count; j++)
                    {
                        if (polylines_2[i][ki].DistanceTo(pts[j]) < 0.001)
                        {
                            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i + nH + nV + n1, ki));
                            //constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i + nH, 2 * ki));
                        }
                    }
                }
            }

            // SOLVER
            solver = new DRRelax(elements, constraints);
            solver.Eclim = Eclim;
            Ec.Clear();
            Ec.Add(solver.Ect);
        }
        private void Init2()
        {
            // reset flags, do first iteration, refresh outputs and restart iterative process
            iteration_current = 1;
            flag = LoopFlag.running;

            // ELEMENT H
            for (int i = 0; i < nH; i++)
            {
                elements[i].Xi = elements[i].X(solver.Xt);
            }

            // ELEMENT V
            for (int j = nH; j < nH + nV; j++)
            {
                elements[j].Xi = elements[j].X(solver.Xt);
            }

            // ELEMENT D1
            for (int j = nH+nV; j < nH + nV+n1; j++)
            {
                elements[j].Xi = elements[j].X(solver.Xt);
            }

            // ELEMENT D2
            for (int j = nH + nV + n1; j < nH + nV + n1 + n2; j++)
            {
                elements[j].Xi = elements[j].X(solver.Xt);
            }
            

            constraints.Clear();
            // LINKS
            int[,] link = new int[4, 2];
            for (int i = 0; i < nH; i++)
            {
                for (int ki = 0; ki < polylines_H[i].Count; ki++)
                {
                    link[0, 0] = i;
                    link[0, 1] = ki;
                    for (int j = 0; j < nV; j++)
                    {
                        for (int kj = 0; kj < polylines_V[j].Count; kj++)
                        {
                            if (polylines_H[i][ki].DistanceTo(polylines_V[j][kj]) < 0.001)
                            {
                                link[1, 0] = j + nH;
                                link[1, 1] = kj;
                                dij = 0;
                                for (int u1 = 0; u1 < n1; u1++)
                                {
                                    for (int k1 = 0; k1 < polylines_1[u1].Count; k1++)
                                    {
                                        if (polylines_H[i][ki].DistanceTo(polylines_1[u1][k1]) < 0.001)
                                        {
                                            dij = 1;
                                            for (int u2 = 0; u2 < n2; u2++)
                                            {
                                                for (int k2 = 0; k2 < polylines_2[u2].Count; k2++)
                                                {
                                                    if (polylines_H[i][ki].DistanceTo(polylines_2[u2][k2]) < 0.001)
                                                    {
                                                        link[2, 0] = u1 + nH + nV;
                                                        link[2, 1] = k1;
                                                        link[3, 0] = u2 + nH + nV + n1;
                                                        link[3, 1] = k2;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                int[] ei = new int[4] { link[0, 0], link[1, 0], link[2, 0], link[3, 0] };
                                int[] nj = new int[4] { link[0, 1], link[1, 1], link[2, 1], link[3, 1] };
                                if (dij == 0)
                                {
                                    constraints.Add(DRLinkConstraint.CreateSingleJointConstraint(i, ki, j + nH, kj));
                                }
                                else
                                {
                                    constraints.Add(DRLinkConstraint.CreateMultipleJointConstraint(ei, nj));
                                }
                            }
                        }

                    }
                }
            }

            // HINGED SUPPORTS
            for (int i = 0; i < nH; i++)
            {
                for (int ki = 0; ki < polylines_H[i].Count; ki++)
                {
                    for (int j = 0; j < pts.Count; j++)
                    {
                        if (polylines_H[i][ki].DistanceTo(pts[j]) < 0.001)
                        {
                            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i, ki));
                            //constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i, 2 * ki));
                        }
                    }
                }
            }
            for (int i = 0; i < nV; i++)
            {
                for (int ki = 0; ki < polylines_V[i].Count; ki++)
                {
                    for (int j = 0; j < pts.Count; j++)
                    {
                        if (polylines_V[i][ki].DistanceTo(pts[j]) < 0.001)
                        {
                            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i + nH, ki));
                            //constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i + nH, 2 * ki));
                        }
                    }
                }
            }
            for (int i = 0; i < n1; i++)
            {
                for (int ki = 0; ki < polylines_1[i].Count; ki++)
                {
                    for (int j = 0; j < pts.Count; j++)
                    {
                        if (polylines_1[i][ki].DistanceTo(pts[j]) < 0.001)
                        {
                            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i + nH + nV, ki));
                            //constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i + nH, 2 * ki));
                        }
                    }
                }
            }
            for (int i = 0; i < n2; i++)
            {
                for (int ki = 0; ki < polylines_2[i].Count; ki++)
                {
                    for (int j = 0; j < pts.Count; j++)
                    {
                        if (polylines_2[i][ki].DistanceTo(pts[j]) < 0.001)
                        {
                            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i + nH + nV + n1, ki));
                            //constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i + nH, 2 * ki));
                        }
                    }
                }
            }

            // SOLVER
            solver = new DRRelax(elements, constraints);
            solver.Eclim = Eclim;
            Ec.Clear();
            Ec.Add(solver.Ect);
        }
        private void Run()
        {
            solver.Run();
            if (solver.Ect < solver.Eclim) { flag = LoopFlag.stoped; }
            Ec.Add(solver.Ect);
        }
        private void OnRefresh()
        {
           
        }
    }
}
