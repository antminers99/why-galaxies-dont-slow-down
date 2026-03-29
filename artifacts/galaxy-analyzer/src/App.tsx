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
