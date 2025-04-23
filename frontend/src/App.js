import logo from './/images/logo.svg';
import './styles/App.css';
import NavBar from './components/NavBar'
import Main from './Main'

function App() {
  return (
      <div className="app">
          <NavBar />
          <Main />
      </div>
  );
}

export default App;
