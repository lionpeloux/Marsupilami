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
    public class GHComponent_2DGridSurfaceAttraction : GH_Component
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
        int nH, nV;
        private List<Curve> crvs_H, crvs_V;
        private Polyline[] polylines_H, polylines_V;
        private Surface srf;
        private List<Point3d> pts;
        private List<double> k, k_cache;
        private double E, S, I, fz;
        bool k_haschanged;


        public DRElement[] elements;
        public List<DRConstraint> constraints;
        public DRRelax solver;
        List<double> Ec;
        public double dt;
        public double Eclim;

        // CONSTRUCTOR
        public GHComponent_2DGridSurfaceAttraction()
            : base("2D grid attraction", "Grid Attraction", "2D grid attraciton by a surface", "Marsupilami", "Example")
        {
            if (DateTime.Now > MarsupilamiInfo.ExpirationDate) { throw new Exception("Enable to launch the gridshell library. Please contact the author."); }

            crvs_H = new List<Curve>();
            crvs_V = new List<Curve>();
            pts = new List<Point3d>();
            k = new List<double>();
            k_cache = new List<double>();

            Ec = new List<double>();
            dt = 1;
            Eclim = 1e-20;

            
        }
        public override Guid ComponentGuid
        {
            get { return new Guid("{002BA6C8-7077-4F63-9C9E-8CFA370344D3}"); }
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
            pManager.AddSurfaceParameter("Surface", "srf", "Attraction surface", GH_ParamAccess.item);
            pManager.AddPointParameter("Attraction Point", "pts", "List of points to attract to the surface", GH_ParamAccess.list);
            pManager.AddNumberParameter("Attractor Stiffness", "k", "List of attraction spring stiffness", GH_ParamAccess.list);
            pManager.AddNumberParameter("Young's Modulus", "E", "Young's modulus of the beam's material (Pa)", GH_ParamAccess.item, 25e9);
            pManager.AddNumberParameter("Cross Section Area", "S", "Area of the beam's cross section (m2)", GH_ParamAccess.item, 4.2e-4);
            pManager.AddNumberParameter("Area Moment of Inertia", "I", "Area moment of inertia of the beam's cross section (m4)", GH_ParamAccess.item, 7.9e-8);
            pManager.AddNumberParameter("Vertical Force", "fz", "Vertical force to apply per node", GH_ParamAccess.item, 0);

            pManager.AddIntegerParameter("iteration_max", "N_max", "Total number of iterations to run", GH_ParamAccess.item, 100000);
            pManager.AddIntegerParameter("iterations_groupby", "N_refresh", "Group of iterations between output refresh", GH_ParamAccess.item, 1000);
            pManager.AddBooleanParameter("reset", "reset", "Reset the engine", GH_ParamAccess.item);

            pManager[0].Optional = false;
            pManager[1].Optional = false;
            pManager[2].Optional = false;
            pManager[3].Optional = false;
            pManager[4].Optional = false;
            pManager[5].Optional = false;
            pManager[6].Optional = false;
            pManager[7].Optional = true;
            pManager[8].Optional = false;
            pManager[9].Optional = false;
            pManager[10].Optional = false;
            pManager[11].Optional = false;
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
            k.Clear();
          
            if (!DA.GetData(2, ref srf)) { return; }
            if (!DA.GetDataList(3, pts)) { return; }
            if (!DA.GetDataList(4, k)) { return; }
            if (!DA.GetData(5, ref E)) { return; }
            if (!DA.GetData(6, ref S)) { return; }
            if (!DA.GetData(7, ref I)) { return; }
            if (!DA.GetData(8, ref fz)) { return; }

            if (!DA.GetData(9, ref iteration_max)) { return; }
            if (!DA.GetData(10, ref iteration_groupby)) { return; }
            if (!DA.GetData(11, ref loop_reset)) { return; }

            
            // DETECT CHANGES IN k
            if (k.Count != k_cache.Count)
            {
                k_haschanged = true;
            }
            else
            {
                k_haschanged = false;
                for (int i = 0; i < k.Count; i++)
                {
                    if (k[i] != k_cache[i])
                    {
                        k_haschanged = true;
                        break;
                    }
                }
            }
                

            if (solver == null || loop_reset != loop_reset_cache) // Premier Calcul
            {
                // EXTRACT POLYLINES FROM INPUT CURVES 
                crvs_H.Clear();
                crvs_V.Clear();
                if (!DA.GetDataList(0, crvs_H)) { return; }
                if (!DA.GetDataList(1, crvs_V)) { return; }

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

                // INIT
                loop_reset_cache = loop_reset;
                Init();

                // FIRST ITERATION
                ghDocument.SolutionEnd += documentSolutionEnd;
                Run();
            }
            else if (k_haschanged)
            {
                Rhino.RhinoApp.WriteLine("k has changed");
                k_cache.Clear(); 
                k_haschanged = false;
                k_cache.AddRange(k);

                // INIT
                loop_reset_cache = loop_reset;
                Init2();

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

            elements = new DRElement[nH + nV];
            constraints = new List<DRConstraint>();
            
            // ELEMENT H
            for (int i = 0; i < nH; i++)
            {
                int n = polylines_H[i].Count;
                var Xi = new Vector3d[n];
                var L0 = new double[n-1];
                var Fext = new Vector3d[n];
                
                Xi[0] = (Vector3d)polylines_H[i][0];
                Fext[0] = new Vector3d(0,0,fz);
                for (int k = 1; k < n; k++)
			    {
                    Xi[k] = (Vector3d)polylines_H[i][k];
                    Fext[k] = new Vector3d(0,0,fz);
                    L0[k-1] = (Xi[k]- Xi[k-1]).Length;
			    }

                elements[i] = DRBeamElement.CreateSimpleBeamElement(i, Xi, L0, Fext, dt);
                elements[i].E = E;
                elements[i].S = S;
                elements[i].I = I;
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
                elements[j].E = E;
                elements[j].S = S;
                elements[j].I = I;
            }


            // LINKS
            for (int i = 0; i < nH; i++)
            {
                for (int j = 0; j < nV; j++)
                {
                    for (int ki = 0; ki < polylines_H[i].Count; ki++)
                    {
                        for (int kj = 0; kj < polylines_V[j].Count; kj++)
                        {
                            if (polylines_H[i][ki].DistanceTo(polylines_V[j][kj]) < 0.001)
                            {
                                constraints.Add(DRLinkConstraint.CreateSingleJointConstraint(i, ki, j+nH, kj));
                            }
                        }
                        
                    }
                }
            }

            // ATTRACTORS
            for (int i = 0; i < nH; i++)
            {
                for (int ki = 0; ki < polylines_H[i].Count; ki++)
                {
                    for (int j = 0; j < pts.Count; j++)
			        {
                        if (polylines_H[i][ki].DistanceTo(pts[j]) < 0.001)
                        {
                            constraints.Add(DRMechanicalConstraint.CreateSurfaceAttractorConstraint(i,ki,srf,k[j]));
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

            // LINKS
            constraints.Clear();
            for (int i = 0; i < nH; i++)
            {
                for (int j = 0; j < nV; j++)
                {
                    for (int ki = 0; ki < polylines_H[i].Count; ki++)
                    {
                        for (int kj = 0; kj < polylines_V[j].Count; kj++)
                        {
                            if (polylines_H[i][ki].DistanceTo(polylines_V[j][kj]) < 0.001)
                            {
                                constraints.Add(DRLinkConstraint.CreateSingleJointConstraint(i, ki, j + nH, kj));
                            }
                        }

                    }
                }
            }

            // ATTRACTORS
            for (int i = 0; i < nH; i++)
            {
                for (int ki = 0; ki < polylines_H[i].Count; ki++)
                {
                    for (int j = 0; j < pts.Count; j++)
                    {
                        if (polylines_H[i][ki].DistanceTo(pts[j]) < 0.001)
                        {
                            constraints.Add(DRMechanicalConstraint.CreateSurfaceAttractorConstraint(i, ki, srf, k[j]));
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
