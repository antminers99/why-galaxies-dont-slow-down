import { Switch, Route, Router as WouterRouter } from "wouter";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { Toaster } from "@/components/ui/toaster";
import { TooltipProvider } from "@/components/ui/tooltip";
import NotFound from "@/pages/not-found";
import Dashboard from "@/pages/dashboard";
import UploadPage from "@/pages/upload";
import AnalysisPage from "@/pages/analysis";
import ModelsPage from "@/pages/models";
import ResearchPage from "@/pages/research";
import CorrelationsPage from "@/pages/correlations";
import TheoryPage from "@/pages/theory";
import StressTestPage from "@/pages/stress-test";
import RARAnalysisPage from "@/pages/rar-analysis";
import ReplicationPage from "@/pages/replication";
import DarkMatterFractionPage from "@/pages/dark-matter-fraction";
import BaryonHaloCouplingPage from "@/pages/baryon-halo-coupling";
import EvidencePage from "@/pages/evidence";
import DefensePage from "@/pages/defense";
import ModelProposalPage from "@/pages/model-proposal";
import ConclusionsPage from "@/pages/conclusions";
import EquationPage from "@/pages/equation";
import RedshiftLabPage from "@/pages/redshift-lab";
import ClusterTestPage from "@/pages/cluster-test";
import BreakTestPage from "@/pages/break-test";
import PipelinePage from "@/pages/pipeline";
import WhyA0Page from "@/pages/why-a0";
import { GalaxyProvider } from "@/hooks/use-galaxy";

const queryClient = new QueryClient();

function Router() {
  return (
    <Switch>
      <Route path="/" component={Dashboard} />
      <Route path="/upload" component={UploadPage} />
      <Route path="/analysis" component={AnalysisPage} />
      <Route path="/models" component={ModelsPage} />
      <Route path="/research" component={ResearchPage} />
      <Route path="/correlations" component={CorrelationsPage} />
      <Route path="/theory" component={TheoryPage} />
      <Route path="/stress-test" component={StressTestPage} />
      <Route path="/rar-analysis" component={RARAnalysisPage} />
      <Route path="/evidence" component={EvidencePage} />
      <Route path="/replication" component={ReplicationPage} />
      <Route path="/dark-matter-fraction" component={DarkMatterFractionPage} />
      <Route path="/baryon-halo-coupling" component={BaryonHaloCouplingPage} />
      <Route path="/defense" component={DefensePage} />
      <Route path="/model" component={ModelProposalPage} />
      <Route path="/conclusions" component={ConclusionsPage} />
      <Route path="/equation" component={EquationPage} />
      <Route path="/redshift-lab" component={RedshiftLabPage} />
      <Route path="/cluster-test" component={ClusterTestPage} />
      <Route path="/break-test" component={BreakTestPage} />
      <Route path="/pipeline" component={PipelinePage} />
      <Route path="/why-a0" component={WhyA0Page} />
      <Route component={NotFound} />
    </Switch>
  );
}

function App() {
  return (
    <QueryClientProvider client={queryClient}>
      <GalaxyProvider>
        <TooltipProvider>
          <WouterRouter base={import.meta.env.BASE_URL.replace(/\/$/, "")}>
            <Router />
          </WouterRouter>
          <Toaster />
        </TooltipProvider>
      </GalaxyProvider>
    </QueryClientProvider>
  );
}

export default App;
