using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Marsupilami.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace Marsupilami.GHComponent
{
    public class GHComponent_CableNet : GH_Component
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
        private List<Point3d> pts;
        private GH_Structure<GH_Vector> Fext_tree, Fext_tree_cache;

        private double ES, ES_cache;
        private bool EC_haschanged, Fext_haschanged;
        private bool is_firstrun;

        public DRElement[] elements;
        public List<DRConstraint> constraints;
        public DRRelax solver;
        List<double> Ec;
        public double dt;
        public double Eclim;

        // CONSTRUCTOR
        public GHComponent_CableNet()
            : base("Cable Net", "Cable Net", "Cable Net", "Marsupilami", "Example")
        {
            if (DateTime.Now > MarsupilamiInfo.ExpirationDate) { throw new Exception("Enable to launch the gridshell library. Please contact the author."); }
            Ec = new List<double>();
            dt = 1;
            Eclim = 1e-16;

            crvs_H = new List<Curve>();
            crvs_V = new List<Curve>();
            pts = new List<Point3d>();

            is_firstrun = true;
        }
        public override Guid ComponentGuid
        {
            get { return new Guid("{1d798c9d-a726-45ad-9860-f5394d395967}"); }
        }
        public override GH_Exposure Exposure
        {
            get { return GH_Exposure.primary; }
        }

        // PARAMETERS
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("H", "H", "", GH_ParamAccess.list);
            pManager.AddCurveParameter("V", "V", "", GH_ParamAccess.list);

            pManager.AddPointParameter("support points", "pts", "Defines the supporting nodes", GH_ParamAccess.list);
            pManager.AddNumberParameter("Axial Stiffness (ES)", "ES", "Cross section axial stifness", GH_ParamAccess.item, 4e-4 * 210e9);
            pManager.AddVectorParameter("Fext", "Fext", "external forces on H nodes (as a tree)", GH_ParamAccess.tree);

            pManager.AddIntegerParameter("iteration_max", "N_max", "Total number of iterations to run", GH_ParamAccess.item, 100000);
            pManager.AddIntegerParameter("iterations_groupby", "N_refresh", "Group of iterations between output refresh", GH_ParamAccess.item, 1000);
            pManager.AddBooleanParameter("reset", "reset", "Reset the engine", GH_ParamAccess.item);

            pManager[0].Optional = false;
            pManager[1].Optional = false;
            pManager[2].Optional = false;
            pManager[3].Optional = false;
            pManager[4].Optional = false;
            pManager[5].Optional = false;
            pManager[6].Optional = true;
            pManager[7].Optional = false;
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
            pManager.Register_IntegerParam("index", "index", "");
        }

        // SOLVER
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            ghDocument = this.OnPingDocument();

            pts.Clear();

           
            if (!DA.GetData(3, ref ES)) { return; }
            if (!DA.GetDataTree(4, out Fext_tree)) { return; }
    
            if (!DA.GetData(5, ref iteration_max)) { return; }
            if (!DA.GetData(6, ref iteration_groupby)) { return; }
            if (!DA.GetData(7, ref loop_reset)) { return; }


            if (!is_firstrun)
            {
                // DETECT CHANGES
                if (ES != ES_cache)
                {
                    EC_haschanged = true;
                    ES_cache = ES;
                }
                Fext_haschanged = TreeHasChanged(Fext_tree, Fext_tree_cache);
                Fext_tree_cache = Fext_tree.Duplicate();
            }

            if (solver == null || loop_reset != loop_reset_cache) // Premier Calcul
            {
                is_firstrun = false;
                Fext_tree_cache = Fext_tree.Duplicate();

                // EXTRACT POLYLINES FROM INPUT CURVES 
                crvs_H.Clear();
                crvs_V.Clear();
                if (!DA.GetDataList(0, crvs_H)) { return; }
                if (!DA.GetDataList(1, crvs_V)) { return; }
                if (!DA.GetDataList(2, pts)) { return; }

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
            else if (EC_haschanged)
            {
                Rhino.RhinoApp.WriteLine("EC has changed");
                EC_haschanged = false;
                
                // INIT
                loop_reset_cache = loop_reset;
                Init_ES();

                // FIRST ITERATION
                ghDocument.SolutionEnd += documentSolutionEnd;
                Run();
            }
            else if (Fext_haschanged)
            {
                Rhino.RhinoApp.WriteLine("Fext has changed");
                Fext_haschanged = false;
                
                // INIT
                loop_reset_cache = loop_reset;
                Init_Fext();

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
                var L0 = new double[n - 1];
          
                var Fext = new Vector3d[n];
                for (int k = 0; k < Fext_tree.Branches[i].Count; k++)
                {
                    Fext[k] = Fext_tree.Branches[i][k].Value;
                }
                
                Xi[0] = (Vector3d)polylines_H[i][0];
                
                for (int k = 1; k < n; k++)
                {
                    Xi[k] = (Vector3d)polylines_H[i][k];
                    L0[k - 1] = (Xi[k] - Xi[k - 1]).Length;
                }

                elements[i] = DRCableElement.CreateSimpleCableElement(i, Xi, L0, Fext, dt);
                elements[i].E = 1;
                elements[i].S = ES;
            }

            // ELEMENT V
            for (int j = nH; j < nH + nV; j++)
            {
                int n = polylines_V[j - nH].Count;
                var Xi = new Vector3d[n];
                var L0 = new double[n - 1];
                var Fext = new Vector3d[n];

                Xi[0] = (Vector3d)polylines_V[j - nH][0];

                for (int k = 1; k < n; k++)
                {
                    Xi[k] = (Vector3d)polylines_V[j - nH][k];
                    L0[k - 1] = (Xi[k] - Xi[k - 1]).Length;
                }

                elements[j] = DRCableElement.CreateSimpleCableElement(j, Xi, L0, Fext, dt);
                elements[j].E = 1;
                elements[j].S = ES;
            }


            // JOINTS
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

            // BOUNDARY CONDITIONS
            for (int i = 0; i < nH; i++)
            {
                for (int ki = 0; ki < polylines_H[i].Count; ki++)
                {
                    for (int j = 0; j < pts.Count; j++)
                    {
                        if (polylines_H[i][ki].DistanceTo(pts[j]) < 0.001)
                        {
                            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(i, ki));
                        }
                    }
                }
            }
            for (int j = 0; j < nV; j++)
            {
                for (int kj = 0; kj < polylines_V[j].Count; kj++)
                {
                    for (int i = 0; i < pts.Count; i++)
                    {
                        if (polylines_V[j][kj].DistanceTo(pts[i]) < 0.001)
                        {
                            constraints.Add(DRKinematicConstraint.CreateSphericalConstraint(j + nH, kj));
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
        private void Init_ES()
        {
            // reset flags, do first iteration, refresh outputs and restart iterative process
            iteration_current = 1;
            flag = LoopFlag.running;

            // ELEMENT H
            for (int i = 0; i < nH; i++)
            {
                elements[i].Xi = elements[i].X(solver.Xt);
                elements[i].S = ES;
            }

            // ELEMENT V
            for (int j = nH; j < nH + nV; j++)
            {
                elements[j].Xi = elements[j].X(solver.Xt);
                elements[j].S = ES;
            }
            
            // SOLVER
            solver = new DRRelax(elements, constraints);
            solver.Eclim = Eclim;
            Ec.Clear();
            Ec.Add(solver.Ect);
        }
        private void Init_Fext()
        {
            // reset flags, do first iteration, refresh outputs and restart iterative process
            iteration_current = 1;
            flag = LoopFlag.running;

            // ELEMENT H
            for (int i = 0; i < nH; i++)
            {
                elements[i].Xi = elements[i].X(solver.Xt);

                for (int k = 0; k < Fext_tree.Branches[i].Count; k++)
                {
                    elements[i].Fext[k] = Fext_tree.Branches[i][k].Value;
                }
            }

            // ELEMENT V
            for (int j = nH; j < nH + nV; j++)
            {
                elements[j].Xi = elements[j].X(solver.Xt);
                elements[j].S = ES;
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

        private bool TreeHasChanged(GH_Structure<GH_Vector> tree, GH_Structure<GH_Vector> tree_cache)
        {
            for (int b = 0; b < tree_cache.Branches.Count; b++)
            {
                for (int i = 0; i < tree_cache.Branches[b].Count; i++)
                {
                    if ((tree_cache.Branches[b][i].Value - tree.Branches[b][i].Value).Length > 0.1)
                    {
                        var v = (tree_cache.Branches[b][i].Value - tree.Branches[b][i].Value).Length;
                        return true;
                    }
                }
            }
            return false;
        }
    }
}
